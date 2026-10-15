//! Test helpers shared by the level-1 routines: random operands, tolerant
//! comparison, and running a check against every core's kernels.

const std = @import("std");
const sysinfo = @import("builtin");
const zigstd = @import("../../../zigstd.zig");
const cpu = @import("../cpu/__root__.zig");
const registry = @import("../asm_registry.zig");
const Vtable = @import("../Vtable.zig");
const blas1 = @import("__root__.zig");
const portable = @import("../portable.zig");
const types = @import("../types.zig");

/// The extended float the `q`/`x` portable kernels are tested at: `XDouble`
/// where it is one, else `f80`, whose kernels are generic Zig ones
/// (`types.kernelCc`) while the vtable's `q`/`x` slots are `f64`.
pub const Ext = if (types.isExtended(types.XDouble)) types.XDouble else f80;

/// `XDouble` and `ComplexXD` where they are types of their own, for the type
/// lists of tests that go through the vtable; empty where they alias `f64`.
pub const x_real = if (types.XDouble == f64) .{} else .{types.XDouble};
pub const x_complex = if (types.XDouble == f64) .{} else .{types.ComplexXD};

/// Lengths that straddle every unroll and block size of the kernels. Not 0:
/// the interface returns before calling a kernel with `n <= 0`, and kernels
/// rely on it (`x86_64/copy_sse2.S` copies an element before testing `n`).
pub const sizes = [_]usize{ 1, 2, 3, 4, 7, 8, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 257 };

/// Element strides: contiguous, and strided in either operand.
pub const strides = [_][2]usize{ .{ 1, 1 }, .{ 2, 1 }, .{ 1, 3 }, .{ 2, 2 } };

/// The longest operand `sizes` and `strides` need, in elements.
pub const capacity = 257 * 3;

pub const isComplex = types.isComplex;
pub const FloatOf = types.Float;

/// Values in [-2, 2), from `seed`; past `f64`, with every significand bit
/// random.
pub fn fill(comptime T: type, buf: []T, seed: u64) void {
    var prng = std.Random.DefaultPrng.init(seed);
    const r = prng.random();
    for (buf) |*v| {
        if (comptime isComplex(T)) {
            v.* = T.init(random(FloatOf(T), r), random(FloatOf(T), r));
        } else {
            v.* = random(T, r);
        }
    }
}

fn random(comptime F: type, r: std.Random) F {
    if (comptime !types.isExtended(F)) return @floatCast(r.float(f64) * 4 - 2);
    // [1, 2) with every fraction bit random, built from bits: converting a
    // u128 is a libcall the self-hosted x86_64 backend leaves on the x87 stack
    const U = @Int(.unsigned, @bitSizeOf(F));
    const frac: U = (@as(U, 1) << std.math.floatFractionalBits(F)) - 1;
    const v: F = @bitCast(@as(U, @bitCast(@as(F, 1))) | (@as(U, @truncate(r.int(u128))) & frac));
    return v * 4 - 6;
}

/// The tolerance for a result of `terms` rounded operations on values of at
/// most 4 in magnitude: kernels sum in any order, the reference in one.
pub fn tolerance(comptime F: type, terms: usize) F {
    // by width: a 64-bit `c_longdouble` is its own type but gets f64's allowance
    const per: F = switch (@typeInfo(F).float.bits) {
        32 => 1e-6,
        64 => 1e-14,
        else => 8 * std.math.floatEps(F),
    };
    return per * 16 * @as(F, @floatFromInt(@max(terms, 1)));
}

pub fn expectClose(comptime T: type, want: T, got: T, tol: FloatOf(T)) !void {
    if (comptime isComplex(T)) {
        try close(want.re, got.re, tol);
        try close(want.im, got.im, tol);
    } else {
        try close(want, got, tol);
    }
}

// Not `std.testing`: its failure message blocks forever on a freestanding rig.
fn close(want: anytype, got: @TypeOf(want), tol: @TypeOf(want)) !void {
    if (std.math.approxEqAbs(@TypeOf(want), want, got, tol)) return;
    say("got {e}, want {e} within {e}: ", .{ got, want, tol });
    return error.TestExpectedApproxEqAbs;
}

/// `expectClose` of `got` against a reference summed in `Wide(F)`; past
/// `f64` compared there, as narrowing `f128` to `f80` is a libcall the
/// self-hosted x86_64 backend leaves on the x87 stack.
pub fn expectCloseRef(comptime F: type, want: Wide(F), got: F, tol: Wide(F)) !void {
    if (comptime types.isExtended(F)) return close(want, @as(Wide(F), got), tol);
    return close(@as(F, @floatCast(want)), got, @as(F, @floatCast(tol)));
}

pub fn expectCloseSlices(comptime T: type, want: []const T, got: []const T, tol: FloatOf(T)) !void {
    try zigstd.testExpectEqual(want.len, got.len);
    for (want, got, 0..) |w, g, i| expectClose(T, w, g, tol) catch |e| {
        say("element {d}: ", .{i});
        return e;
    };
}

// ---- a wide reference, for any float type ----------------------------------

/// `T` at the precision of the `f128` reference.
pub fn Wide(comptime T: type) type {
    return if (isComplex(T)) std.math.Complex(f128) else f128;
}

pub fn wide(comptime T: type, v: T) Wide(T) {
    return if (comptime isComplex(T)) .init(v.re, v.im) else v;
}

/// `|re| + |im|` of a `Wide` value, `|v|` for a real one.
pub fn abs1(v: anytype) f128 {
    return if (comptime isComplex(@TypeOf(v))) @abs(v.re) + @abs(v.im) else @abs(v);
}

pub fn mulW(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.mul(b) else a * b;
}

pub fn addW(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.add(b) else a + b;
}

pub fn conjW(a: anytype) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.conjugate() else a;
}

/// `got` against the reference `want`, within `ulps` units of `T`'s
/// precision of `scale`: a backward-error bound, `scale` being the sum of the
/// magnitudes that went into `want`.
pub fn expectNear(comptime T: type, want: Wide(T), got: T, ulps: f128, scale: f128) !void {
    const F = FloatOf(T);
    const tol = ulps * std.math.floatEps(F) * scale + std.math.floatMin(F);
    try expectClose(Wide(T), want, wide(T, got), tol);
}

/// The float types beyond the vtable's to test the interface at: `f80` too
/// where it may cross the C ABI.
pub const any_floats = if (sysinfo.cpu.arch.isX86()) .{ f16, f80, f128 } else .{ f16, f128 };

/// Says which case of a check failed: on stderr, or a qemu rig's console.
pub fn say(comptime f: []const u8, args: anytype) void {
    const root = @import("root");
    if (comptime sysinfo.os.tag != .freestanding) {
        std.debug.print(f, args);
    } else if (comptime @hasDecl(root, "testWrite")) {
        var buf: [256]u8 = undefined;
        root.testWrite(std.fmt.bufPrint(&buf, f, args) catch return);
    }
}

/// Names a failed case of `routine` before passing its error on.
pub fn case(e: anyerror, routine: []const u8, comptime T: type, n: usize, inc: [2]usize) anyerror {
    say("{s} {s} n={d} inc={d},{d}: ", .{ routine, @typeName(T), n, inc[0], inc[1] });
    return e;
}

/// Runs `check` on a vtable holding each linked core's kernels over the
/// portable ones, for every core this machine runs. A failing core does not
/// stop the rest, so one run names them all.
pub fn eachCore(comptime check: fn (*const Vtable) anyerror!void) !void {
    var failed: ?anyerror = null;
    inline for (@typeInfo(cpu.CoreType).@"enum".field_names) |name| {
        const core = @field(cpu.CoreType, name);
        if (registry.isLinked(core) and (!@hasDecl(cpu.arch, "canRun") or cpu.arch.canRun(core))) {
            var vt: Vtable = .{};
            portable.install(&vt);
            cpu.kernel.install(&vt, core);
            portable.finish(&vt);
            check(&vt) catch |e| {
                say("core {s}\n", .{name});
                failed = e;
            };
        }
    }
    if (failed) |e| return e;
}

test "tolerance goes by width, so a 64-bit c_longdouble gets f64's" {
    const L = c_longdouble;
    const want: L = if (@typeInfo(L).float.bits == 64) 1e-14 * 16 else 8 * std.math.floatEps(L) * 16;
    try zigstd.testExpectEqual(want, tolerance(L, 1));
    try zigstd.testExpectEqual(8 * std.math.floatEps(f80) * 16, tolerance(f80, 1));
}
