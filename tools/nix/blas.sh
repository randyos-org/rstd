#!/usr/bin/env bash
# nixpkgs' prebuilt OpenBLAS and BLIS for blas_bench's library sections
# (build/sshtest.zig), from the nixpkgs of flake.lock; substituted from the
# binary cache, never built (--max-jobs 0). Run on the Nix machine by
# tools/ssh_test.zig.
#
#   blas.sh <system> <nixpkgs rev> [<bench nix dir>]
#
# With the bench's nix directory (tools/blasbench/nix), also the reference
# binaries (tools/blasbench/readme.md): Reference-LAPACK (lapack-reference),
# AOCL-libFLAME (amd-libflame, x86_64-linux only; skipped when the cache
# lacks it), and the Eigen shim, which this machine builds (eigen.nix).
#
# Linux: writes a newc cpio to stdout of the libraries' runtime closure (each
# store path's lib/*.so*), /flyoff/lib{openblas,blis,lapack,eigen,flame}.so,
# and glibc's loader where Zig's executables look for it. Darwin: links them
# under /tmp/rstd-blas/<system> and prints their paths.
set -euo pipefail
sys=$1
np=github:NixOS/nixpkgs/$2
bench=${3:-}
export PATH="/nix/var/nix/profiles/default/bin:$PATH"

# openblasCompat: 32-bit integers (plain openblas is ILP64 on x86_64-linux);
# nixpkgs builds BLIS for x86_64-linux only
pkgs=(openblasCompat)
if [ "$sys" = x86_64-linux ]; then pkgs+=(blis); fi
if [ -n "$bench" ]; then pkgs+=(lapack-reference); fi
# the store path of package $1; $2 an out-link (a GC root), else none
real() {
  local link=(--no-link)
  if [ $# -gt 1 ]; then link=(--out-link "$2"); fi
  nix build --max-jobs 0 --print-out-paths "${link[@]}" "$np#legacyPackages.$sys.$1.out"
}
# the Eigen shim for $sys, built here by nixpkgs' Zig; $1 an out-link, else none.
# The CPU is the flyoff host's for each system (readme.md).
eigen() {
  local link=(--no-link) here target mcpu
  if [ $# -gt 0 ]; then link=(--out-link "$1"); fi
  here=$(nix eval --impure --raw --expr builtins.currentSystem)
  case $sys in
  x86_64-linux) target=x86_64-linux-gnu mcpu=znver3 ;;
  aarch64-linux) target=aarch64-linux-gnu mcpu=apple_m1 ;;
  aarch64-darwin) target=aarch64-macos mcpu=apple_m1 ;;
  *)
    echo "blas.sh: no Eigen target for $sys" >&2
    return 2
    ;;
  esac
  nix build --impure --print-out-paths "${link[@]}" --expr \
    "(builtins.getFlake \"$np\").legacyPackages.$here.callPackage (/. + \"$bench/eigen.nix\") { target = \"$target\"; mcpu = \"$mcpu\"; }"
}

case $sys in
*-darwin)
  out=/tmp/rstd-blas/$sys
  mkdir -p "$out"
  for p in "${pkgs[@]}"; do
    s=$(real "$p" "$out/$p")
    echo "$p: $out/$p -> $s"
  done
  if [ -n "$bench" ]; then
    s=$(eigen "$out/eigen")
    echo "eigen: $out/eigen -> $s"
  fi
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
t=$(mktemp -d)
trap 'rm -rf "$t"' EXIT
mkdir -p "$t/flyoff" "$t/${ld%/*}"
# each package's root, and its library under /flyoff
link() { # <store path> <library in its lib/> <name under /flyoff>
  roots+=("$1")
  ln -s "$1/lib/$2" "$t/flyoff/$3"
}
# (assigned first: `set -e` sees a failed substitution only there)
p=$(real openblasCompat)
link "$p" libopenblas.so libopenblas.so
if [ "$sys" = x86_64-linux ]; then
  p=$(real blis)
  link "$p" libblis.so libblis.so
fi
if [ -n "$bench" ]; then
  p=$(real lapack-reference)
  link "$p" liblapack.so liblapack.so
  p=$(eigen)
  link "$p" libeigen_shim.so libeigen.so
  if [ "$sys" = x86_64-linux ]; then
    if flame=$(real amd-libflame); then
      link "$flame" libflame.so libflame.so
    else
      echo "blas.sh: amd-libflame is not in the binary cache; its column is skipped" >&2
    fi
  fi
fi
glibc=$(real glibc)
roots+=("$glibc")
ln -s "$glibc/lib/${ld##*/}" "$t/$ld"

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
