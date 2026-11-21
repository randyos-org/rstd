#!/usr/bin/env bash
# kernel.sh <system> <nixpkgs rev> <x86_64-darwin rev>: nixpkgs' prebuilt Linux kernel for <system> (cache only, --max-jobs 0) on
# stdout, for the Linux rigs and the qemu-bench-linux steps (build/tools.zig's nixKernel, issue #656). Run on the Nix machine, as
# blas.sh is; an x86_64-darwin one evaluates the release that still supports it (flake.nix).
set -euo pipefail
export PATH="/nix/var/nix/profiles/default/bin:$PATH"
rev=$2
if [ "$(nix eval --impure --raw --expr builtins.currentSystem)" = x86_64-darwin ]; then rev=$3; fi
# the cache holds only these two; the other systems' kernels would build from source
case $1 in
x86_64-linux) image=bzImage ;;
aarch64-linux) image=Image ;;
*)
  echo "kernel.sh: no cached nixpkgs kernel for $1" >&2
  exit 2
  ;;
esac
out=$(nix build --max-jobs 0 --no-link --print-out-paths "github:NixOS/nixpkgs/$rev#legacyPackages.$1.linux.out")
cat "$out/$image"
