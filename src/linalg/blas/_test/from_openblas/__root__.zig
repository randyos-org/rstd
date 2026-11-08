//! OpenBLAS's own suites (`test/[sdcz]blat[123].f` with their `.dat` files,
//! `ctest/`, `utest/`, `test/compare_*.c`), OpenBLAS v0.3.34-317 (a69062524);
//! tables and data transcribed mechanically: re-derive them, don't edit (`design.md`'s Verification).

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../types.zig");
const check = @import("../check.zig");
const blas = @import("../../__root__.zig");
const Vtable = @import("../../Vtable.zig");
const util = @import("../../util.zig");
const blat1 = @import("blat1.zig");
const blat2 = @import("blat2.zig");
const blat3 = @import("blat3.zig");
const utest = @import("utest.zig");
const utest_ext = @import("utest_ext.zig");
const compare = @import("compare.zig");

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

/// A check's test ratio, and whether the result is less than half accurate.
pub const Result = struct { err: f64, fatal: bool };

/// A leading dimension one past `ld`, as the Fortran sweeps take it, up to
/// `max`; null past it.
pub fn bump(ld: usize, max: usize) ?usize {
    const l = if (ld < max) ld + 1 else ld;
    return if (l > max) null else l;
}

/// The `Trans` of a Fortran letter: `N`, `T`, else `C`.
pub fn transOf(c: u8) blas.Trans {
    return switch (c) {
        'N' => .n,
        'T' => .t,
        else => .c,
    };
}

pub fn sideOf(c: u8) blas.Side {
    return if (c == 'L') .l else .r;
}

pub fn uploOf(c: u8) blas.Uplo {
    return if (c == 'U') .u else .l;
}

pub fn diagOf(c: u8) blas.Diag {
    return if (c == 'U') .u else .n;
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

pub const val = check.val;

pub const zero = util.zero;

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
    return if (comptime isComplex(@TypeOf(v))) util.abs(v.re) + util.abs(v.im) else util.abs(v);
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

/// `suite` on `vt`'s kernels, the sweeps on the rigs' table. The suites go
/// in pieces, so each test fits a rig image (`build/qemu/machines.zig`'s
/// `sparc64-debug`).
fn OnCore(comptime suite: anytype) fn (*const Vtable) anyerror!void {
    return struct {
        fn run(vt: *const Vtable) anyerror!void {
            const saved = blas.testUse(vt);
            defer blas.testRestore(saved);
            try suite();
        }
    }.run;
}

fn level1Utest() anyerror!void {
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
    try utest_ext.level3All();
    try compare.all(true);
}

test "blat1 on every linked core this machine runs" {
    try check.eachCore(OnCore(blat1.all));
}

test "utest and utest_ext's level-1 cases on every linked core this machine runs" {
    try check.eachCore(OnCore(level1Utest));
}

test "utest_ext's level-2 cases on every linked core this machine runs" {
    try check.eachCore(OnCore(utest_ext.level2All));
}

test "blat2 at f32 on every linked core this machine runs" {
    try check.eachCore(OnCore(Level2(f32).run));
}

test "blat2 at f64 on every linked core this machine runs" {
    try check.eachCore(OnCore(Level2(f64).run));
}

test "blat2 at complex f32 on every linked core this machine runs" {
    try check.eachCore(OnCore(Level2(std.math.Complex(f32)).run));
}

test "blat2 at complex f64 on every linked core this machine runs" {
    try check.eachCore(OnCore(Level2(std.math.Complex(f64)).run));
}

test "blat3, utest's LAPACK cases, utest_ext's level-3 cases and the comparisons on every linked core this machine runs" {
    try check.eachCore(OnCore(level3));
}

test {
    _ = blat1;
    _ = blat2;
    _ = blat3;
    _ = utest;
    _ = utest_ext;
    _ = compare;
}
