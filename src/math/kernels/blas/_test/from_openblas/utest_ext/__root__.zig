//! OpenBLAS's extension unit tests, `utest/test_extensions/test_*.c` and
//! `test_bgemm.c`, OpenBLAS v0.3.34-317; each file's tables at its end are the
//! sources' cases. Seeded operands, not `rand()`; norms count equal infinities and NaNs as agreeing.

const vector = @import("vector.zig");
const level1 = @import("level1.zig");
const level2 = @import("level2.zig");
const geadd = @import("geadd.zig");
const matcopy = @import("matcopy.zig");
const gemm = @import("gemm.zig");
const bgemm = @import("bgemm.zig");

/// The vector cases, for `check.eachCore`; with `level2All` and `level3All`,
/// every case here. Apart, so each fits a rig image.
pub fn all() !void {
    try vector.reductions();
    try level1.axpby();
    try level1.axpyc();
    try level1.complexRot();
    try level1.complexScal();
    try level1.complexRotg();
    try level1.rotmg();
}

/// The level-2 cases, for `check.eachCore`.
pub fn level2All() !void {
    try level2.gbmv();
    try level2.gemvOS();
    try level2.gemvTCUD();
    try level2.sbmvSpmv();
    try level2.trmvTrsv();
}

/// The matrix extensions and `gemm`, for `check.eachCore`; they run with the
/// level-3 suite, whose drivers they reach.
pub fn level3All() !void {
    try geadd.geadd();
    try matcopy.matcopy();
    try gemm.complexGemm();
    try gemm.gemmt();
    try bgemm.bgemm();
}

test {
    _ = vector;
    _ = level1;
    _ = level2;
    _ = geadd;
    _ = matcopy;
    _ = gemm;
    _ = bgemm;
}
