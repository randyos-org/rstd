//! OpenBLAS's own test suites, ported to rstd's interfaces: the reference
//! BLAS test programs (`test/[sdcz]blat[123].f` and their `.dat` inputs),
//! their CBLAS versions (`ctest/`) and the unit tests (`utest/`), all from
//! OpenBLAS v0.3.34-317 (commit a69062524). Each file names its sources. The
//! sweep tables, data and expected values are transcribed mechanically; if
//! they go stale, re-derive them from those files rather than edit them.
//!
//! The heavy sweeps run their `.dat` tables natively; on the rigs and WASI
//! (`reduced`) they take the smaller table each program's header documents.
//! One more test runs every suite, on that table, on each linked core the
//! machine runs (`check.eachCore`).
//!
//! Not translatable, so not ported:
//! - `xerbla` itself, and the error exits for a bad letter argument: rstd
//!   takes enums, so only the letters a routine rejects (`O` to `D` but in
//!   complex `gemv` and `gbmv`, `V` and `M` but in the Hermitian level-2
//!   routines) stand in for them;
//! - checks that scalar arguments come back unchanged (Zig passes them by
//!   value), and the snapshot, summary and stop-on-failure options;
//! - CBLAS: the row-major layout and the argument numbers shifted by the
//!   order argument; the column-major ctest cases are the Fortran ones;
//! - `utest/test_fork.c`, `test_post_fork*.c` (processes and threads), and the
//!   `ARMV8` `sdot`/`ddot` cases, which preload a vector register by inline
//!   assembly before the call.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../types.zig");
const check = @import("../blas1/check.zig");
const blas = @import("../../__root__.zig");
const Vtable = @import("../../Vtable.zig");

pub const Float = types.Float;
pub const isComplex = types.isComplex;
pub const say = check.say;
/// The `Io` the suites pass the routines.
pub const io = check.io;

/// The value of an argument enum (`blas.Trans` and the like) a Fortran letter
/// names, either case; null for a letter it has none for.
pub fn letter(comptime E: type, c: u8) ?E {
    return std.meta.stringToEnum(E, &.{std.ascii.toLower(c)});
}

/// On a rig: freestanding (32 MiB heap, 1 MiB stack, emulated), WASI, or the
/// Linux rigs, whose init sets `RSTD_TEST_PROGRESS`.
pub fn reduced() bool {
    return switch (sysinfo.target.os.tag) {
        .linux => std.testing.environ.getPosix("RSTD_TEST_PROGRESS") != null,
        .freestanding, .other, .uefi, .wasi => true,
        else => false,
    };
}

/// The precision a check computes its expected values in: `f64` for single
/// precision, the working precision otherwise, as the Fortran does.
pub fn Ref(comptime T: type) type {
    const R = if (Float(T) == f32) f64 else Float(T);
    return if (isComplex(T)) std.math.Complex(R) else R;
}

pub fn val(comptime T: type, re: f64, im: f64) T {
    return if (comptime isComplex(T)) .init(@floatCast(re), @floatCast(im)) else @floatCast(re);
}

pub fn zero(comptime T: type) T {
    return val(T, 0, 0);
}

pub fn isZero(v: anytype) bool {
    return if (comptime isComplex(@TypeOf(v))) v.re == 0 and v.im == 0 else v == 0;
}

pub fn eql(a: anytype, b: @TypeOf(a)) bool {
    return if (comptime isComplex(@TypeOf(a))) a.re == b.re and a.im == b.im else a == b;
}

pub fn add(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.add(b) else a + b;
}

pub fn sub(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.sub(b) else a - b;
}

pub fn mul(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.mul(b) else a * b;
}

pub fn conj(a: anytype) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.conjugate() else a;
}

/// `|re| + |im|`, the Fortran's `ABS1`; `|v|` for a real `v`.
pub fn abs1(v: anytype) Float(@TypeOf(v)) {
    return if (comptime isComplex(@TypeOf(v))) @abs(v.re) + @abs(v.im) else @abs(v);
}

/// The modulus, the Fortran's `ABS`.
pub fn absv(v: anytype) Float(@TypeOf(v)) {
    return if (comptime isComplex(@TypeOf(v))) std.math.hypot(v.re, v.im) else @abs(v);
}

/// `v` in the type `U` of the same kind.
pub fn cast(comptime U: type, v: anytype) U {
    if (comptime isComplex(U)) return .init(@floatCast(v.re), @floatCast(v.im));
    return @floatCast(v);
}

/// The Fortran's `SBEG`/`CBEG`: values in [-0.5, 0.5) from a fixed sequence
/// that repeats with a period of 50 (`RESET` is a fresh `Gen`).
pub const Gen = struct {
    i: i32 = 7,
    j: i32 = 7,
    ic: i32 = 0,

    pub fn next(g: *Gen, comptime T: type) T {
        g.ic += 1;
        while (true) {
            g.i = @mod(g.i * 891, 1000);
            g.j = @mod(g.j * 457, 1000);
            if (g.ic < 5) break;
            g.ic = 0;
        }
        const F = Float(T);
        const re = @as(F, @floatFromInt(g.i - 500)) / 1001;
        if (comptime !isComplex(T)) return re;
        return .init(re, @as(F, @floatFromInt(g.j - 500)) / 1001);
    }
};

/// The Fortran's `LSE`: whether two arrays are identical.
pub fn lse(comptime T: type, a: []const T, b: []const T) bool {
    for (a, b) |p, q| if (!eql(p, q)) return false;
    return true;
}

/// Fails, naming `what`, unless `ok`. Not `std.testing`: its failure message
/// blocks forever on a freestanding rig.
pub fn expect(ok: bool, what: []const u8) !void {
    if (ok) return;
    say("{s}\n", .{what});
    return error.TestUnexpectedResult;
}

/// Fails unless `got` is the error `want`; `what` names the case.
pub fn expectError(want: anyerror, got: anytype, what: []const u8) !void {
    if (got) |_| {
        say("{s}: no error, want {s}\n", .{ what, @errorName(want) });
        return error.TestExpectedError;
    } else |e| if (e != want) {
        say("{s}: {s}, want {s}\n", .{ what, @errorName(e), @errorName(want) });
        return error.TestUnexpectedError;
    }
}

const blat1 = @import("blat1.zig");
const blat2 = @import("blat2.zig");
const blat3 = @import("blat3.zig");
const utest = @import("utest.zig");
const utest_ext = @import("utest_ext.zig");

/// `suite` on `vt`'s kernels, the sweeps on the rigs' table. The suites go
/// in pieces, so each test fits a rig image (`build/qemutest.zig`'s
/// `sparc64-debug`).
fn OnCore(comptime suite: fn () anyerror!void) fn (*const Vtable) anyerror!void {
    return struct {
        fn run(vt: *const Vtable) anyerror!void {
            const saved = blas.testUse(vt);
            defer blas.testRestore(saved);
            try suite();
        }
    }.run;
}

fn level1() anyerror!void {
    try blat1.all();
    try utest.all();
    try utest_ext.all();
}

fn Level2(comptime T: type) type {
    return struct {
        fn run() anyerror!void {
            try blat2.allOf(T, &blat2.rig_n);
        }
    };
}

fn level3() anyerror!void {
    try blat3.all(&blat2.rig_n);
    try utest.lapackAll();
}

test "the level-1 suites on every linked core this machine runs" {
    try check.eachCore(OnCore(level1));
}

test "the level-2 suite at f32 on every linked core this machine runs" {
    try check.eachCore(OnCore(Level2(f32).run));
}

test "the level-2 suite at f64 on every linked core this machine runs" {
    try check.eachCore(OnCore(Level2(f64).run));
}

test "the level-2 suite at complex f32 on every linked core this machine runs" {
    try check.eachCore(OnCore(Level2(std.math.Complex(f32)).run));
}

test "the level-2 suite at complex f64 on every linked core this machine runs" {
    try check.eachCore(OnCore(Level2(std.math.Complex(f64)).run));
}

test "the level-3 suite on every linked core this machine runs" {
    try check.eachCore(OnCore(level3));
}

test {
    _ = blat1;
    _ = blat2;
    _ = blat3;
    _ = utest;
    _ = utest_ext;
}
