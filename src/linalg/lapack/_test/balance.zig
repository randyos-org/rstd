//! Tests of `../balance.zig`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const ops = @import("../ops.zig");
const blas_calls = @import("../blas_calls.zig");
const balancing = @import("../balance.zig");

const Allocator = std.mem.Allocator;
const Float = types.Float;
const isComplex = types.isComplex;
const Side = types.Side;
const BalanceJob = balancing.BalanceJob;
const Block = balancing.Block;
const back = balancing.back;
const balance = balancing.balance;
const gebak = balancing.gebak;
const gebal = balancing.gebal;

const testing = std.testing;
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");

const test_types = lcheck.test_types;

/// `TESTING/dbal.in` of Reference-LAPACK 3.12.0 (OpenBLAS v0.3.34-317's
/// `lapack-netlib/`), verbatim: `dchkbl.f`'s cases for `xGEBAL` with job
/// `B`, each `n`, `A` by rows, the expected `ilo`, `ihi`, balanced `A` and
/// `scale`. Copy it again rather than editing it if it goes stale.
const dbal_in =
    \\DBL:  Tests DGEBAL
    \\  5
    \\  0.1000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.2000D+01  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.3000D+01  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.4000D+01  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.5000D+01
    \\
    \\   1   1
    \\  0.1000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.2000D+01  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.3000D+01  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.4000D+01  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.5000D+01
    \\
    \\  0.1000D+01  0.2000D+01  0.3000D+01  0.4000D+01  0.5000D+01
    \\
    \\  5
    \\  0.1000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.1000D+01  0.2000D+01  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.1000D+01  0.2000D+01  0.3000D+01  0.0000D+00  0.0000D+00
    \\  0.1000D+01  0.2000D+01  0.3000D+01  0.4000D+01  0.0000D+00
    \\  0.1000D+01  0.2000D+01  0.3000D+01  0.4000D+01  0.5000D+01
    \\
    \\   1   1
    \\  0.5000D+01  0.4000D+01  0.3000D+01  0.2000D+01  0.1000D+01
    \\  0.0000D+00  0.4000D+01  0.3000D+01  0.2000D+01  0.1000D+01
    \\  0.0000D+00  0.0000D+00  0.3000D+01  0.2000D+01  0.1000D+01
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.2000D+01  0.1000D+01
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.1000D+01
    \\
    \\  0.1000D+01  0.2000D+01  0.3000D+01  0.2000D+01  0.1000D+01
    \\
    \\  5
    \\  0.1000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.1000D+01  0.1000D+01  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.1000D+01  0.1000D+01  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.1000D+01  0.1000D+01  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.1000D+01  0.1000D+01
    \\
    \\   1   1
    \\  0.1000D+01  0.1000D+01  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.1000D+01  0.1000D+01  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.1000D+01  0.1000D+01  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.1000D+01  0.1000D+01
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.1000D+01
    \\
    \\  0.1000D+01  0.2000D+01  0.3000D+01  0.2000D+01  0.1000D+01
    \\
    \\  4
    \\  0.0000D+00  0.2000D+01  0.1000D+00  0.0000D+00
    \\  0.2000D+01  0.0000D+00  0.0000D+00  0.1000D+00
    \\  0.1000D+03  0.0000D+00  0.0000D+00  0.2000D+01
    \\  0.0000D+00  0.1000D+03  0.2000D+01  0.0000D+00
    \\
    \\   1   4
    \\  0.0000D-03  2.0000D+00  3.2000D+00  0.0000D-03
    \\  2.0000D+00  0.0000D-03  0.0000D-03  3.2000D+00
    \\  3.1250D+00  0.0000D-03  0.0000D-03  2.0000D+00
    \\  0.0000D-03  3.1250D+00  2.0000D+00  0.0000D-03
    \\
    \\  62.5000D-03 62.5000D-03  2.0000D+00  2.0000D+00
    \\
    \\  6
    \\  0.2000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.1024D+04
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.1280D+03
    \\  0.0000D+00  0.2000D+01  0.3000D+04  0.0000D+00  0.0000D+00  0.2000D+01
    \\  0.1280D+03  0.4000D+01  0.4000D-02  0.5000D+01  0.6000D+03  0.8000D+01
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.2000D-02  0.2000D+01
    \\  0.8000D+01  0.8192D+04  0.0000D+00  0.0000D+00  0.0000D+00  0.2000D+01
    \\
    \\   4   6
    \\  0.5000D+01  0.4000D-02  0.6000D+03  0.1024D+04  0.5000D+00  0.8000D+01
    \\  0.0000D+00  0.3000D+04  0.0000D+00  0.0000D+00  0.2500D+00  0.2000D+01
    \\  0.0000D+00  0.0000D+00  0.2000D-02  0.0000D+00  0.0000D+00  0.2000D+01
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.2000D+01  0.0000D+00  0.1280D+03
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.1024D+04
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.6400D+02  0.1024D+04  0.2000D+01
    \\
    \\  0.4000D+01  0.3000D+01  0.5000D+01  0.8000D+01  0.1250D+00  0.1000D+01
    \\
    \\  5
    \\  0.1000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.8000D+01
    \\  0.0000D+00  0.2000D+01  0.8192D+04  0.2000D+01  0.4000D+01
    \\  0.2500D-03  0.1250D-03  0.4000D+01  0.0000D+00  0.6400D+02
    \\  0.0000D+00  0.2000D+01  0.1024D+04  0.4000D+01  0.8000D+01
    \\  0.0000D+00  0.8192D+04  0.0000D+00  0.0000D+00  0.8000D+01
    \\
    \\   1   5
    \\  1.0000D+00     0.0000D-03     0.0000D-03     0.0000D-03     2.0000D+00
    \\  0.0000D-03     2.0000D+00     1.0240D+03    16.0000D+00    16.0000D+00
    \\  3.2000D-02     1.0000D-03     4.0000D+00     0.0000D-03     2.0480D+03
    \\  0.0000D-03   250.0000D-03    16.0000D+00     4.0000D+00     4.0000D+00
    \\  0.0000D-03     2.0480D+03     0.0000D-03     0.0000D-03     8.0000D+00
    \\
    \\  8.0000D+00  500.0000D-03  62.5000D-03  4.0000D+00  2.0000D+00
    \\
    \\  4
    \\  0.1000D+01  0.1000D+07  0.1000D+07  0.1000D+07
    \\ -0.2000D+07  0.3000D+01  0.2000D-05  0.3000D-05
    \\ -0.3000D+07  0.0000D+00  0.1000D-05  0.2000D+01
    \\  0.1000D+07  0.0000D+00  0.3000D-05  0.4000D+07
    \\
    \\   1   4
    \\  1.0000D+00     1.0000D+06     2.0000D+06     1.0000D+06
    \\ -2.0000D+06     3.0000D+00     4.0000D-06     3.0000D-06
    \\ -1.5000D+06     0.0000D-03     1.0000D-06     1.0000D+00
    \\  1.0000D+06     0.0000D-03     6.0000D-06     4.0000D+06
    \\
    \\  1.0000D+00  1.0000D+00 2.0000D+00  1.0000D+00
    \\
    \\   4
    \\  0.1000D+01  0.1000D+05  0.1000D+05  0.1000D+05
    \\ -0.2000D+05  0.3000D+01  0.2000D-02  0.3000D-02
    \\  0.0000D+00  0.2000D+01  0.0000D+00 -0.3000D+05
    \\  0.0000D+00  0.0000D+00  0.1000D+05  0.0000D+00
    \\
    \\   1   4
    \\  1.0000D+00    10.0000D+03    10.0000D+03     5.0000D+03
    \\ -20.0000D+03     3.0000D+00     2.0000D-03     1.5000D-03
    \\  0.0000D-03     2.0000D+00     0.0000D-03   -15.0000D+03
    \\  0.0000D-03     0.0000D-03    20.0000D+03     0.0000D-03
    \\
    \\   1.0000D+00     1.0000D+00     1.0000D+00   500.0000D-03
    \\
    \\  5
    \\  0.1000D+01  0.5120D+03  0.4096D+04  3.2768D+04  2.62144D+05
    \\  0.8000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.8000D+01  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.8000D+01  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.8000D+01  0.0000D+00
    \\
    \\   1   5
    \\  1.0000D+00    32.0000D+00   32.0000D+00  32.0000D+000 32.0000D+00
    \\  128.0000D+00  0.0000D-03    0.0000D-03   0.0000D-003  0.0000D-03
    \\  0.0000D-03    64.0000D+00   0.0000D-03   0.0000D-003  0.0000D-03
    \\  0.0000D-03    0.0000D-03    64.0000D+00  0.0000D-003  0.0000D-03
    \\  0.0000D-03    0.0000D-03    0.0000D-03   64.0000D+000 0.0000D-03
    \\
    \\  256.0000D+00  16.0000D+00  2.0000D+00  250.0000D-03  31.2500D-03
    \\
    \\  6
    \\  0.1000D+01  0.1000D+01  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.1000D+01  0.1000D+01  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.1000D+01  0.0000D+00  0.0000D+00
    \\  0.1000D+01  0.1000D+01  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.1000D+01  0.1000D+01  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01
    \\
    \\   2   5
    \\  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.1000D+01
    \\
    \\  0.3000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.4000D+01
    \\
    \\  7
    \\  0.6000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.1000D+01 0.0000D+00
    \\  0.0000D+00  0.4000D+01  0.0000D+00  0.2500D-03  0.1250D-01  0.2000D-01 0.1250D+00
    \\  0.1000D+01  0.1280D+03  0.6400D+02  0.0000D+00  0.0000D+00 -0.2000D+01 0.1600D+02
    \\  0.0000D+00  1.6384D+04  0.0000D+00  0.1000D+01 -0.4000D+03  0.2560D+03 -0.4000D+04
    \\ -0.2000D+01 -0.2560D+03  0.0000D+00  0.1250D-01  0.2000D+01  0.2000D+01 0.3200D+02
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00 0.0000D+00
    \\  0.0000D+00  0.8000D+01  0.0000D+00  0.4000D-02  0.1250D+00 -0.2000D+00 0.3000D+01
    \\
    \\  2   5
    \\  6.4000D+01   1.0000D+00   5.00000D-01   0.0000D+00   0.0000D+00   1.0000D+00  -2.0000D+00
    \\  0.0000D+00   4.0000D+00   5.00000D-01   1.0240D+00   8.0000D-01   0.0000D+00   2.5600D+00
    \\  0.0000D+00   2.0000D+00   3.00000D+00   4.0960D+00   2.0000D+00   0.0000D+00  -6.4000D+00
    \\  0.0000D+00   4.0000D+00  -3.90625D+00   1.0000D+00  -6.2500D+00   0.0000D+00   8.0000D+00
    \\  0.0000D+00  -4.0000D+00   2.00000D+00   8.0000D-01   2.0000D+00  -4.0000D+00   4.0000D+00
    \\  0.0000D+00   0.0000D+00   0.00000D+00   0.0000D+00   0.0000D+00   6.0000D+00   1.0000D+00
    \\  0.0000D+00   0.0000D+00   0.00000D+00   0.0000D+00   0.0000D+00   0.0000D+00   0.0000D+00
    \\
    \\  3.0000D+00  7.812500D-03  3.1250D-02  3.2000D+01  5.0000D-01  1.0000D+00 6.0000D+00
    \\
    \\  5
    \\  0.1000D+04  0.2000D+01  0.3000D+01  0.4000D+01  0.5000D+06
    \\  0.9000D+01  0.0000D+00  0.2000D-03  0.1000D+01  0.3000D+01
    \\  0.0000D+00 -0.3000D+03  0.2000D+01  0.1000D+01  0.1000D+01
    \\  0.9000D+01  0.2000D-02  0.1000D+01  0.1000D+01 -0.1000D+04
    \\  0.6000D+01  0.2000D+03  0.1000D+01  0.6000D+03  0.3000D+01
    \\
    \\  1   5
    \\  1.0000D+03     3.1250D-02     3.7500D-01     3.1250D-02     1.95312500D+03
    \\  5.7600D+02     0.0000D+00     1.6000D-03     5.0000D-01     7.50000000D-01
    \\  0.0000D+00    -3.7500D+01     2.0000D+00     6.2500D-02     3.12500000D-02
    \\  1.1520D+03     4.0000D-03     1.6000D+01     1.0000D+00    -5.00000000D+02
    \\  1.5360D+03     8.0000D+02     3.2000D+01     1.2000D+03     3.00000000D+00
    \\
    \\  3.2000D+01  5.0000D-01   4.0000D+00   2.5000D-01  1.2500D-01
    \\
    \\  6
    \\  1.0000D+00  1.0000D+120 0.0000D+00  0.0000D+00  0.0000D+00 0.0000D+00
    \\  1.0000D-120 1.0000D+00  1.0000D+120 0.0000D+00  0.0000D+00 0.0000D+00
    \\  0.0000D+00  1.0000D-120 1.0000D+00  1.0000D+120 0.0000D+00 0.0000D+00
    \\  0.0000D+00  0.0000D+00  1.0000D-120 1.0000D+00  1.0000D+120 0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  1.0000D-120 1.0000D+00 1.0000D+120
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  1.0000D-120 1.0000D+00
    \\
    \\   1   6
    \\0.10000000000000000000D+01 0.63448545932891229313D+04 0.00000000000000000000D+00 0.00000000000000000000D+00 0.00000000000000000000D+00 0.00000000000000000000D+00
    \\0.15760802478557791348D-03 0.10000000000000000000D+01 0.63448545932891229313D+04 0.00000000000000000000D+00 0.00000000000000000000D+00 0.00000000000000000000D+00
    \\0.00000000000000000000D+00 0.15760802478557791348D-03 0.10000000000000000000D+01 0.31724272966445614657D+04 0.00000000000000000000D+00 0.00000000000000000000D+00
    \\0.00000000000000000000D+00 0.00000000000000000000D+00 0.31521604957115582695D-03 0.10000000000000000000D+01 0.15862136483222807328D+04 0.00000000000000000000D+00
    \\0.00000000000000000000D+00 0.00000000000000000000D+00 0.00000000000000000000D+00 0.63043209914231165391D-03 0.10000000000000000000D+01 0.79310682416114036641D+03
    \\0.00000000000000000000D+00 0.00000000000000000000D+00 0.00000000000000000000D+00 0.00000000000000000000D+00 0.12608641982846233078D-02 0.10000000000000000000D+01
    \\
    \\  2.494800386918399765D+291  1.582914569427869018D+175  1.004336277661868922D+59  3.186183822264904554D-58  5.053968264940243633D-175  0.40083367200179455560D-291
    \\
    \\0
;

/// The next number of `dbal_in`, Fortran's `D` exponent read as `e`.
fn nextNumber(it: *std.mem.TokenIterator(u8, .any)) f64 {
    const tok = it.next().?;
    var buf: [64]u8 = undefined;
    const s = buf[0..tok.len];
    for (tok, s) |c, *d| d.* = if (c == 'D' or c == 'd') 'e' else c;
    return std.fmt.parseFloat(f64, s) catch unreachable;
}

/// `gebak` against `gebal`'s `A` and `B = D⁻¹ Pᵀ A P D`: `A X = X B` with
/// `X = P D` (right), `Yᵀ A = B Yᵀ` with `Y = P D⁻¹` (left), within
/// `4 n eps ‖A‖ ‖X‖` per entry.
fn expectBack(comptime T: type, al: Allocator, job: BalanceJob, n: usize, a: []const T, b: []const T, blk: Block, scale: []const Float(T)) !void {
    const F = Float(T);
    const x = try al.alloc(T, n * n);
    defer al.free(x);
    var amax: F = 0;
    for (a, b) |p, q| amax = @max(amax, ops.abs(T, p), ops.abs(T, q));
    for ([_]Side{ .r, .l }) |side| {
        @memset(x, ops.zero(T));
        for (0..n) |i| x[i * (n + 1)] = ops.one(T);
        try gebak(check.io, T, job, side, n, blk.lo, blk.hi, scale, n, x, n);
        var xmax: F = 0;
        for (x) |v| xmax = @max(xmax, ops.abs(T, v));
        const tol = 4 * @as(F, @floatFromInt(n)) * std.math.floatEps(F) * amax * xmax;
        for (0..n) |i| for (0..n) |j| {
            var l = ops.zero(T);
            var r = ops.zero(T);
            for (0..n) |k| {
                if (side == .r) {
                    l = ops.add(T, l, ops.mul(T, a[i + k * n], x[k + j * n]));
                    r = ops.add(T, r, ops.mul(T, x[i + k * n], b[k + j * n]));
                } else {
                    l = ops.add(T, l, ops.mul(T, x[k + i * n], a[k + j * n]));
                    r = ops.add(T, r, ops.mul(T, b[i + k * n], x[j + k * n]));
                }
            }
            const d = ops.abs(T, ops.sub(T, l, r));
            if (!(d <= tol)) {
                check.say("gebak {s} job {s} side {s} n {d}: ({d}, {d}) off by {e}, tolerance {e}\n", .{ @typeName(T), @tagName(job), @tagName(side), n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(tol)) });
                return error.TestUnexpectedResult;
            }
        };
    }
}

test "lapack gebal: the reference's dbal.in cases, and gebak takes vectors back" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        var it = std.mem.tokenizeAny(u8, dbal_in, " \n\r");
        // the title line's three words
        for (0..3) |_| _ = it.next();
        var case: usize = 1;
        while (true) : (case += 1) {
            const n: usize = @intFromFloat(nextNumber(&it));
            if (n == 0) break;
            var a: [7 * 7]f64 = undefined;
            var want: [7 * 7]f64 = undefined;
            var wscale: [7]f64 = undefined;
            for (0..n) |i| for (0..n) |j| {
                a[i + j * n] = nextNumber(&it);
            };
            // the file's 1-based `ilo`, `ihi` and swapped rows
            const wlo: usize = @as(usize, @intFromFloat(nextNumber(&it))) - 1;
            const whi: usize = @intFromFloat(nextNumber(&it));
            for (0..n) |i| for (0..n) |j| {
                want[i + j * n] = nextNumber(&it);
            };
            for (0..n) |i| wscale[i] = nextNumber(&it);
            // case 13's 1e120 and its scalings are past f32's range
            var fits = true;
            for (a[0 .. n * n]) |v| if (v != 0 and (@abs(v) > std.math.floatMax(F) or @abs(v) < std.math.floatMin(F))) {
                fits = false;
            };
            for (wscale[0..n]) |v| if (@abs(v) > std.math.floatMax(F) or @abs(v) < std.math.floatMin(F)) {
                fits = false;
            };
            if (!fits) continue;
            var got: [7 * 7]T = undefined;
            for (0..n * n) |i| got[i] = lcheck.fromF64(T, a[i]);
            const a0 = got;
            var scale: [7]F = undefined;
            const blk = try gebal(check.io, T, .b, n, &got, n, &scale);
            try testing.expectEqual(wlo, blk.lo);
            try testing.expectEqual(whi, blk.hi);
            // A and scale checked only by expectBack for case 5, whose swap
            // order is not the one 3.12.0's `dgebal.f` loops (ported line for
            // line) give, and for case 13 in a type past f64's range, which
            // scales further
            const exact = case != 5 and !(case == 13 and std.math.floatMax(F) > std.math.floatMax(f64));
            // dchkbl.f's relative difference; the file's digits are exact
            // but for case 13's twenty
            const tol: F = @max(8 * std.math.floatEps(F), 1e-15);
            if (exact) for (0..n * n) |i| {
                const g = got[i];
                const w: F = @floatCast(want[i]);
                const d = if (comptime isComplex(T)) @abs(g.re - w) + @abs(g.im) else @abs(g - w);
                if (!(d <= tol * @max(ops.abs(T, g), @abs(w), std.math.floatMin(F)))) {
                    check.say("gebal {s} dbal.in case {d}: A[{d}] {e} against {e}\n", .{ @typeName(T), case, i, @as(f64, @floatCast(ops.re(T, g))), want[i] });
                    return error.TestUnexpectedResult;
                }
            };
            if (exact) for (0..n) |i| {
                const w: F = @floatCast(wscale[i] - @as(f64, if (i < wlo or i >= whi) 1 else 0));
                if (!(@abs(scale[i] - w) <= tol * @max(scale[i], w))) {
                    check.say("gebal {s} dbal.in case {d}: scale[{d}] {e} against {e}\n", .{ @typeName(T), case, i, @as(f64, @floatCast(scale[i])), wscale[i] });
                    return error.TestUnexpectedResult;
                }
            };
            try expectBack(T, al, .b, n, a0[0 .. n * n], got[0 .. n * n], blk, scale[0..n]);
        }
    }
}

test "lapack gebal: every job on random and structured matrices, with gebak" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 5, 12 }) |n| for ([_]BalanceJob{ .n, .p, .s, .b }) |job| for (0..3) |kind| {
            const a = try al.alloc(T, n * n);
            defer al.free(a);
            check.fill(T, a, n * 31 + kind);
            // badly scaled rows and columns
            for (0..n) |i| for (0..n) |j| {
                const e: F = @floatFromInt(@as(isize, @intCast(i)) - @as(isize, @intCast(j)));
                a[i + j * n] = ops.mulReal(T, a[i + j * n], @as(F, @floatCast(std.math.pow(f64, 4, @floatCast(e)))));
            };
            // isolated eigenvalues: a row, then a column, zero off the diagonal
            if (kind >= 1 and n > 2) for (0..n) |j| if (j != 1) {
                a[1 + j * n] = ops.zero(T);
            };
            if (kind == 2 and n > 3) for (0..n) |i| if (i != n - 2) {
                a[i + (n - 2) * n] = ops.zero(T);
            };
            const b = try al.alloc(T, n * n);
            defer al.free(b);
            @memcpy(b, a);
            const scale = try al.alloc(F, n);
            defer al.free(scale);
            const blk = try gebal(check.io, T, job, n, b, n, scale);
            try testing.expect(blk.lo <= blk.hi and blk.hi <= n);
            if (job == .n) {
                try testing.expect(blk.lo == 0 and blk.hi == n);
                try testing.expectEqualSlices(u8, std.mem.sliceAsBytes(a), std.mem.sliceAsBytes(b));
            }
            if (job == .s) try testing.expect(blk.lo == 0 and blk.hi == n);
            // below the diagonal outside lo..hi: zero
            for (0..n) |j| for (j + 1..n) |i| {
                if (j < blk.lo or i >= blk.hi) try testing.expect(ops.isZero(T, b[i + j * n]));
            };
            // the scalings are powers of 2
            if (job == .s or job == .b) for (blk.lo..blk.hi) |i| {
                try testing.expectEqual(@as(F, 0.5), std.math.frexp(scale[i]).significand);
            };
            try expectBack(T, al, job, n, a, b, blk, scale);
        };
    }
}

test "lapack gebal: n = 0 and 1, NaN, bad arguments; gebak's" {
    inline for (test_types) |T| {
        const F = Float(T);
        var scale = [_]F{ 7, 7 };
        var a = [_]T{ ops.scalar(T, 2), ops.scalar(T, 1), ops.scalar(T, 1), ops.scalar(T, 3) };
        var blk = try gebal(check.io, T, .b, 0, &a, 1, &scale);
        try testing.expect(blk.lo == 0 and blk.hi == 0 and scale[0] == 7);
        blk = try gebal(check.io, T, .b, 1, &a, 1, &scale);
        try testing.expect(blk.lo == 0 and blk.hi == 1 and scale[0] == 1);
        a[1] = lcheck.poisoned(T);
        try testing.expectError(error.NanInA, gebal(check.io, T, .s, 2, &a, 2, &scale));
    }
    var a = [_]f64{ 1, 2, 3, 4 };
    var s = [_]f64{ 1, 1 };
    try testing.expectError(error.LeadingDimensionTooSmall, gebal(check.io, f64, .b, 2, &a, 1, &s));
    try testing.expectError(error.SliceTooShort, gebal(check.io, f64, .b, 2, a[0..3], 2, &s));
    try testing.expectError(error.SliceTooShort, gebal(check.io, f64, .b, 2, &a, 2, s[0..1]));
    try testing.expectError(error.BlockOutOfRange, gebak(check.io, f64, .b, .r, 1, 1, 0, &s, 1, &a, 1));
    try testing.expectError(error.BlockOutOfRange, gebak(check.io, f64, .b, .l, 2, 0, 3, &s, 1, &a, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, gebak(check.io, f64, .b, .l, 2, 0, 2, &s, 1, &a, 1));
    try testing.expectError(error.SliceTooShort, gebak(check.io, f64, .b, .l, 2, 0, 2, &s, 3, &a, 2));
    try testing.expectError(error.SliceTooShort, gebak(check.io, f64, .b, .l, 2, 0, 2, s[0..1], 1, &a, 2));
    // a swapped row past `V`'s, outside the block
    s = .{ 2, 1 };
    try testing.expectError(error.PivotOutOfRange, gebak(check.io, f64, .p, .l, 2, 1, 2, &s, 1, &a, 2));
    s = .{ 0.5, 1 };
    try testing.expectError(error.PivotOutOfRange, gebak(check.io, f64, .b, .l, 2, 1, 2, &s, 1, &a, 2));
    // job s reads no rows from scale
    try gebak(check.io, f64, .s, .l, 2, 1, 2, &s, 1, &a, 2);
    try gebak(check.io, f64, .b, .l, 0, 0, 0, &s, 1, &a, 1);
}

test "lapack gebal: D⁻¹ Pᵀ A P D is the balanced matrix, isolated eigenvalues outside ilo..ihi" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| {
        const F = Float(T);
        const n = 9;
        const a0 = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, a0);
        check.fill(T, a0.ptr[0 .. n * n], 4);
        // graded rows and columns, and two isolated eigenvalues: row 2 and
        // column 6 zero off the diagonal
        for (0..n) |j| for (0..n) |i| {
            const d = @as(F, @floatCast(std.math.pow(f64, 2, @as(f64, @floatFromInt(i)) * 3 - @as(f64, @floatFromInt(j)) * 3)));
            a0.at(i, j).* = ops.mulReal(T, a0.get(i, j), d);
            if ((i == 2 and j != 2) or (j == 6 and i != 6)) a0.at(i, j).* = ops.zero(T);
        };
        const b = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, b);
        b.copyFrom(a0);
        var scale: [n]F = undefined;
        const blk = try balance(check.io, T, b, &scale);
        try testing.expect(blk.lo >= 1 and blk.hi <= n - 1);
        // X = P D, the identity backed: then A X = X B
        const x = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, x);
        x.identity();
        try back(check.io, T, blk, &scale, false, x);
        const ax = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, ax);
        const xb = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, xb);
        try blas_calls.gemm(al, check.io, T, .n, .n, ops.one(T), a0, x, ops.zero(T), ax);
        try blas_calls.gemm(al, check.io, T, .n, .n, ops.one(T), x, b, ops.zero(T), xb);
        for (0..n) |j| for (0..n) |i| {
            const err = ops.abs(T, ops.sub(T, ax.get(i, j), xb.get(i, j)));
            try testing.expect(err <= 16 * std.math.floatEps(F) * @max(ops.abs(T, ax.get(i, j)), std.math.floatMin(F)));
        };
        // outside the window: triangular
        for (0..n) |j| for (0..n) |i| {
            if (i > j and (j < blk.lo or i >= blk.hi)) try testing.expect(ops.isZero(T, b.get(i, j)));
        };
    }
}
