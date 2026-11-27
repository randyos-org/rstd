#!/usr/bin/env bash
# hintrefs.sh <nixpkgs rev> <tools/hintrefs dir> <font> [axis values, 16.16]: tools/hintrefs/hintrefs.c built against that nixpkgs' FreeType, run on
# <font>, gzipped on stdout (rstd issue #805). Run on a Nix machine by tools/ssh_test.zig; `zig build hint-refs -Dtest-fonts`.
set -euo pipefail
np=github:NixOS/nixpkgs/$1
dir=$2
font=$3
shift 3
export PATH="/nix/var/nix/profiles/default/bin:$PATH"
here=$(nix eval --impure --raw --expr builtins.currentSystem)
exe=$(nix build --impure --no-link --print-out-paths --expr \
  "(builtins.getFlake \"$np\").legacyPackages.$here.callPackage (/. + \"$dir/hintrefs.nix\") { }")
"$exe" "$font" "$@" | gzip -9n
