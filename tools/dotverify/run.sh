#!/bin/sh
# run.sh <OpenBLAS checkout> [zig], from rstd's root: the x86_64 dot drivers (`_verify_dot.zig`) against OpenBLAS's kernel/x86_64/*dot.c,
# bit for bit (T38), per core this CPU runs, fused and not. Files in `$DOTVERIFY_OUT` or a temporary directory; exits 1 on a mismatch.
set -e
# a directory as the compiler names it (MSYS shells: `C:/...`, not `/c/...`)
native() { (cd "$1" && (pwd -W 2>/dev/null || pwd)); }
OB=$(native "${1:?usage: run.sh <OpenBLAS checkout> [zig]}")
ZIG=${2:-zig}
HERE=$(native "$(dirname "$0")")
OUT=${DOTVERIFY_OUT:-$(mktemp -d)}
mkdir -p "$OUT"
OUT=$(native "$OUT")
ARGS=""
for fz in f n; do
  if [ $fz = f ]; then contract=on; else contract=off; fi
  for row in "dsdot sdot.c -DDSDOT" "ddot ddot.c -DDOUBLE" "cdotu cdot.c -DCOMPLEX" \
    "cdotc cdot.c -DCOMPLEX -DCONJ" "zdotu zdot.c -DCOMPLEX -DDOUBLE" "zdotc zdot.c -DCOMPLEX -DDOUBLE -DCONJ"; do
    set -- $row
    slot=$1 src=$2
    shift 2
    # `_verify_dot.zig`'s cores, by the reference's defines
    for core in nehalem sandybridge haswell zen skylakex generic; do
      case $core in
        generic) def="" ;;
        *) def="-D$(echo "$core" | tr a-z A-Z)" ;;
      esac
      name=ref_${slot}_${core}_$fz
      printf '#include "%s/kernel/x86_64/%s"\n' "$OB" "$src" > "$OUT/$name.c"
      ARGS="$ARGS -cflags -DCNAME=$name $def -ffp-contract=$contract -O2 -I $HERE $* -- $OUT/$name.c"
    done
  done
done
"$ZIG" build-exe -lc -O ReleaseFast -mcpu=native $ARGS \
  --dep blaskern_linked --dep blas_options -Mroot=src/linalg/blas/_verify_dot.zig -Mblaskern_linked="$HERE/stub_linked.zig" \
  -Mblas_options="$HERE/stub_options.zig" \
  --name verify_dot --cache-dir "$OUT/cache" -femit-bin="$OUT/verify_dot"
"$OUT/verify_dot"
