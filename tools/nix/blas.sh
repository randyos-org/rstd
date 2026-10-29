#!/usr/bin/env bash
# nixpkgs' prebuilt OpenBLAS and BLIS for blas_bench's library sections
# (build/sshtest.zig), from the nixpkgs of flake.lock; substituted from the
# binary cache, never built (--max-jobs 0). Run on the Nix machine by
# tools/ssh_test.zig.
#
#   blas.sh <system> <nixpkgs rev>
#
# Linux: writes a newc cpio to stdout of the libraries' runtime closure (each
# store path's lib/*.so*), /flyoff/lib{openblas,blis}.so, and glibc's loader
# where Zig's executables look for it. Darwin: links them under
# /tmp/rstd-blas/<system> and prints their paths.
set -euo pipefail
sys=$1
np=github:NixOS/nixpkgs/$2
export PATH="/nix/var/nix/profiles/default/bin:$PATH"

# openblasCompat: 32-bit integers (plain openblas is ILP64 on x86_64-linux);
# nixpkgs builds BLIS for x86_64-linux only
pkgs=(openblasCompat)
if [ "$sys" = x86_64-linux ]; then pkgs+=(blis); fi
# the store path of package $1; $2 an out-link (a GC root), else none
real() {
  local link=(--no-link)
  if [ $# -gt 1 ]; then link=(--out-link "$2"); fi
  nix build --max-jobs 0 --print-out-paths "${link[@]}" "$np#legacyPackages.$sys.$1.out"
}

case $sys in
*-darwin)
  out=/tmp/rstd-blas/$sys
  mkdir -p "$out"
  for p in "${pkgs[@]}"; do
    echo "$p: $out/$p -> $(real "$p" "$out/$p")"
  done
  exit 0
  ;;
esac

case $sys in
x86_64-linux) ld=lib64/ld-linux-x86-64.so.2 ;;
aarch64-linux) ld=lib/ld-linux-aarch64.so.1 ;;
*)
  echo "blas.sh: no loader path for $sys" >&2
  exit 2
  ;;
esac
roots=()
for p in "${pkgs[@]}" glibc; do roots+=("$(real "$p")"); done
glibc=${roots[${#roots[@]} - 1]}

t=$(mktemp -d)
trap 'rm -rf "$t"' EXIT
mkdir -p "$t/flyoff" "$t/${ld%/*}"
ln -s "$glibc/lib/${ld##*/}" "$t/$ld"
ln -s "${roots[0]}/lib/libopenblas.so" "$t/flyoff/libopenblas.so"
if [ ${#roots[@]} -gt 2 ]; then ln -s "${roots[1]}/lib/libblis.so" "$t/flyoff/libblis.so"; fi

# directories before what they hold, names relative to /
{
  printf 'nix\nnix/store\n'
  for s in $(nix-store -qR "${roots[@]}"); do
    printf '%s\n' "${s#/}"
    if [ -d "$s/lib" ]; then
      printf '%s\n' "${s#/}/lib"
      (cd / && find "${s#/}/lib" -maxdepth 1 \( -type f -o -type l \) -name '*.so*')
    fi
  done
} | (cd / && cpio -o -H newc 2>/dev/null)
(cd "$t" && find flyoff "${ld%%/*}" | cpio -o -H newc 2>/dev/null)
