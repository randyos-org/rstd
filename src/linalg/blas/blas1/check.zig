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
const vector = @import("../vector.zig");
const extremum = @import("extremum.zig");

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

/// `strides`, signed, and pairs with exactly one zero, some beside a negative
/// one: element 0 of a zero-stride operand takes every update in turn, which
/// a kernel that loads a block of it before storing gets wrong.
pub const zero_strides = signed: {
    var s: [strides.len][2]isize = undefined;
    for (strides, &s) |p, *q| q.* = .{ p[0], p[1] };
    break :signed s;
} ++ [_][2]isize{ .{ 0, 1 }, .{ 1, 0 }, .{ 0, 2 }, .{ 3, 0 }, .{ -1, 0 }, .{ 0, -1 }, .{ -2, 0 }, .{ 0, -3 } };

/// Keeps `r`'s error in `failed`, so the checks after it still run and one
/// run names every failing kernel.
pub fn keep(failed: *?anyerror, r: anyerror!void) void {
    r catch |e| {
        failed.* = e;
    };
}

/// Where element `i` of a vector of `n` at stride `inc` sits: from the far
/// end for a negative stride, where the interface points the kernel.
pub fn pos(n: usize, inc: isize, i: usize) usize {
    const s: usize = @abs(inc);
    return if (inc < 0) ((n -| 1) - i) * s else i * s;
}

/// `expectCloseSlices` within `tol`, or where `inc` has a zero stride, within
/// a bound for `n` updates in turn of one element, relative past `n`.
pub fn expectStrided(comptime T: type, want: []const T, got: []const T, tol: FloatOf(T), n: usize, inc: [2]isize) !void {
    if (inc[0] != 0 and inc[1] != 0) return expectCloseSlices(T, want, got, tol);
    const F = FloatOf(T);
    const k: F = @floatFromInt(n);
    for (want, got, 0..) |w, g, i| {
        const size: F = if (comptime isComplex(T)) @abs(w.re) + @abs(w.im) else @abs(w);
        expectClose(T, w, g, @max(tol, 16 * k * @max(k, size) * std.math.floatEps(F))) catch |e| {
            say("element {d}: ", .{i});
            return e;
        };
    }
}

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
// NaN matches NaN: an untouched or propagated NaN is what a check expects.
fn close(want: anytype, got: @TypeOf(want), tol: @TypeOf(want)) !void {
    if (std.math.approxEqAbs(@TypeOf(want), want, got, tol)) return;
    if (std.math.isNan(want) and std.math.isNan(got)) return;
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

/// `@sqrt` of a `Wide` value. sparc64's compiler-rt `_Qp_sqrt` takes the root
/// in `f64` (`issues.md`); two Newton steps there restore `f128`'s precision.
pub fn sqrtW(x: f128) f128 {
    var y = @sqrt(x);
    if (comptime sysinfo.cpu.arch != .sparc64) return y;
    if (!(x > 0) or !std.math.isFinite(x)) return y;
    // the guess is taken in f64: keep `x` in its range
    if (x < 0x1p-1000) return sqrtW(x * 0x1p1000) * 0x1p-500;
    if (x > 0x1p1000) return sqrtW(x * 0x1p-1000) * 0x1p500;
    inline for (0..2) |_| y = (y + x / y) / 2;
    return y;
}

/// `|re| + |im|` of a `Wide` value, `|v|` for a real one.
pub fn abs1(v: anytype) f128 {
    return if (comptime isComplex(@TypeOf(v))) @abs(v.re) + @abs(v.im) else @abs(v);
}

/// `abs1` of the finite parts: a tolerance's scale, where the results a
/// special part reaches are NaN or infinite and need none.
pub fn mag(v: anytype) f128 {
    if (comptime isComplex(@TypeOf(v))) return mag(v.re) + mag(v.im);
    return if (std.math.isFinite(v)) @abs(v) else 0;
}

pub fn mulW(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.mul(b) else a * b;
}

pub fn addW(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.add(b) else a + b;
}

pub fn subW(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.sub(b) else a - b;
}

pub fn divW(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) a.div(b) else a / b;
}

/// `f` times `a`, part by part: a real factor, with no cross terms.
pub fn scaleW(f: f128, a: anytype) @TypeOf(a) {
    return if (comptime isComplex(@TypeOf(a))) .init(f * a.re, f * a.im) else f * a;
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

// ---- special values ---------------------------------------------------------

pub const Special = enum { nan, inf, ninf };

/// What a case puts in its random operands.
pub const Pattern = union(enum) {
    /// One special value at element `at`, in its real or imaginary part by
    /// `at`'s parity.
    one: struct { kind: Special, at: usize },
    /// By index mod 8: -0, +0, a denormal, a huge value, +Inf, -Inf, NaN, and
    /// one left random.
    mixed,
    /// Every element denormal.
    denormal,
    /// Every element near the square root of the smallest normal: squares
    /// underflow.
    tiny,
    /// Every element huge, though no product with a value of at most 2, or
    /// sum of 256 of them, overflows; squares do.
    huge,

    pub fn describe(p: Pattern) void {
        switch (p) {
            .one => |o| say("{s} at {d}: ", .{ @tagName(o.kind), o.at }),
            else => say("{s}: ", .{@tagName(p)}),
        }
    }
};

/// The patterns for a vector of `n`: each special value at the first, middle
/// and last element, then the whole-vector ones.
pub fn patterns(n: usize, buf: *[13]Pattern) []const Pattern {
    var len: usize = 0;
    for (std.enums.values(Special)) |kind| for ([_]usize{ 0, n / 2, n - 1 }, 0..) |at, k| {
        if (k > 0 and at == 0 or k == 2 and at == n / 2) continue;
        buf[len] = .{ .one = .{ .kind = kind, .at = at } };
        len += 1;
    };
    for ([_]Pattern{ .mixed, .denormal, .tiny, .huge }) |p| {
        buf[len] = p;
        len += 1;
    }
    return buf[0..len];
}

pub fn special(comptime F: type, kind: Special) F {
    return switch (kind) {
        .nan => std.math.nan(F),
        .inf => std.math.inf(F),
        .ninf => -std.math.inf(F),
    };
}

/// The factors that take a value of at most 2 into each range.
fn denormalScale(comptime F: type) F {
    return comptime std.math.ldexp(@as(F, 1), std.math.floatExponentMin(F) - 4);
}
fn tinyScale(comptime F: type) F {
    return comptime std.math.ldexp(@as(F, 1), @divFloor(std.math.floatExponentMin(F), 2) - 4);
}
fn hugeScale(comptime F: type) F {
    return comptime std.math.ldexp(@as(F, 1), std.math.floatExponentMax(F) - 12);
}

fn scaled(comptime T: type, v: T, f: FloatOf(T)) T {
    return if (comptime isComplex(T)) .init(v.re * f, v.im * f) else v * f;
}

/// `f` in the real part of `v` for an even `i`, else the imaginary part.
fn put(comptime T: type, v: T, i: usize, f: FloatOf(T)) T {
    if (comptime !isComplex(T)) return f;
    return if (i % 2 == 0) .init(f, v.im) else .init(v.re, f);
}

/// Random element `v` at index `i` of an operand, shaped by `p`.
pub fn shape(comptime T: type, v: T, i: usize, p: Pattern) T {
    const F = FloatOf(T);
    return switch (p) {
        .one => |o| if (i == o.at) put(T, v, i, special(F, o.kind)) else v,
        .mixed => switch (i % 8) {
            0, 1 => scaled(T, v, if (i % 8 == 0) -0.0 else 0),
            2 => scaled(T, v, denormalScale(F)),
            3 => scaled(T, v, hugeScale(F)),
            4 => put(T, v, i / 8, std.math.inf(F)),
            5 => put(T, v, i / 8, -std.math.inf(F)),
            6 => put(T, v, i / 8, std.math.nan(F)),
            else => v,
        },
        .denormal => scaled(T, v, denormalScale(F)),
        .tiny => scaled(T, v, tinyScale(F)),
        .huge => scaled(T, v, hugeScale(F)),
    };
}

/// `got` against the reference `want` where special values may appear: NaN
/// for NaN, the infinity of `want`'s sign where it is one or past `T`'s
/// range, else within `ulps` units of `scale` plus `ulps` of `T`'s smallest
/// denormal, so a flushed denormal shows.
pub fn expectWide(comptime T: type, want: Wide(T), got: T, ulps: f128, scale: f128) !void {
    if (comptime isComplex(T)) {
        try component(FloatOf(T), want.re, got.re, ulps, scale);
        try component(FloatOf(T), want.im, got.im, ulps, scale);
    } else try component(T, want, got, ulps, scale);
}

/// `expectWide` of one float with `F`'s range and precision; `got` may be
/// wider (`dsdot`).
pub fn component(comptime F: type, want: f128, got: anytype, ulps: f128, scale: f128) !void {
    if (componentOk(F, want, got, ulps, scale)) return;
    const g: f128 = got;
    say("got {e}, want {e}: ", .{ g, want });
    return error.TestExpectedApproxEqAbs;
}

/// `expectWide` where complex infinities may meet in another order than the
/// reference's: a complex part that should be infinite or NaN may come out
/// either (Inf times a finite complex value is NaN in one order of the
/// multiplications, infinite in another).
pub fn expectWideAnyOrder(comptime T: type, want: Wide(T), got: T, ulps: f128, scale: f128) !void {
    if (comptime isComplex(T)) {
        const lax = !std.math.isFinite(want.re) and !std.math.isFinite(got.re);
        const lay = !std.math.isFinite(want.im) and !std.math.isFinite(got.im);
        if (!lax) try component(FloatOf(T), want.re, got.re, ulps, scale);
        if (!lay) try component(FloatOf(T), want.im, got.im, ulps, scale);
        return;
    }
    return expectWide(T, want, got, ulps, scale);
}

/// Whether `expectWide` would pass, saying nothing.
pub fn wideOk(comptime T: type, want: Wide(T), got: T, ulps: f128, scale: f128) bool {
    if (comptime isComplex(T)) return componentOk(FloatOf(T), want.re, got.re, ulps, scale) and componentOk(FloatOf(T), want.im, got.im, ulps, scale);
    return componentOk(T, want, got, ulps, scale);
}

fn componentOk(comptime F: type, want: f128, got: anytype, ulps: f128, scale: f128) bool {
    const g: f128 = got;
    if (std.math.isNan(want)) return std.math.isNan(g);
    if (@abs(want) > std.math.floatMax(F)) return std.math.isInf(g) and (g > 0) == (want > 0);
    return std.math.isFinite(g) and @abs(g - want) <= ulps * (std.math.floatEps(F) * scale + std.math.floatTrueMin(F));
}

/// Past a vector's last element, NaN that must stay as it was.
pub const tail = 4;

fn Bits(comptime F: type) type {
    return @Int(.unsigned, @bitSizeOf(F));
}

/// The value bits of `a` and `b` agree (not an `f80`'s padding); with
/// `signed_zero` false, zeros of either sign agree too.
pub fn same(comptime T: type, a: T, b: T, signed_zero: bool) bool {
    if (comptime isComplex(T)) return same(FloatOf(T), a.re, b.re, signed_zero) and same(FloatOf(T), a.im, b.im, signed_zero);
    if (!signed_zero and a == 0 and b == 0) return true;
    return @as(Bits(T), @bitCast(a)) == @as(Bits(T), @bitCast(b));
}

/// All ones: what `Vector` and `Matrix` fill unread storage with.
pub fn poisoned(comptime T: type, a: T) bool {
    if (comptime isComplex(T)) return poisoned(FloatOf(T), a.re) and poisoned(FloatOf(T), a.im);
    return @as(Bits(T), @bitCast(a)) == std.math.maxInt(Bits(T));
}

/// A vector operand for edge cases: its `n` elements at stride `inc`, NaN in
/// every gap and for `tail` past the end, and the elements as they were.
pub fn Vector(comptime T: type) type {
    return struct {
        const Self = @This();
        buf: []T,
        was: []T,
        n: usize,
        inc: isize,

        /// Random elements from `seed`, shaped by `p`.
        pub fn init(al: std.mem.Allocator, n: usize, inc: isize, seed: u64, p: ?Pattern) !Self {
            const buf = try al.alloc(T, (n - 1) * @abs(inc) + 1 + tail);
            errdefer al.free(buf);
            var self: Self = .{ .buf = buf, .was = try al.alloc(T, n), .n = n, .inc = inc };
            @memset(std.mem.sliceAsBytes(self.buf), 0xff);
            fill(T, self.was, seed);
            if (p) |q| for (self.was, 0..) |*v, i| {
                v.* = shape(T, v.*, i, q);
            };
            for (0..n) |i| self.buf[self.idx(i)] = self.was[i];
            return self;
        }

        pub fn deinit(self: *Self, al: std.mem.Allocator) void {
            al.free(self.was);
            al.free(self.buf);
        }

        /// Where element `i` lives: from the far end for a negative stride.
        pub fn idx(self: *const Self, i: usize) usize {
            const s: usize = @abs(self.inc);
            return if (self.inc > 0) i * s else (self.n - 1 - i) * s;
        }

        pub fn at(self: *const Self, i: usize) T {
            return self.buf[self.idx(i)];
        }

        /// Element `i` as it was, widened.
        pub fn w(self: *const Self, i: usize) Wide(T) {
            return wide(T, self.was[i]);
        }

        /// Element 0, where a kernel takes the vector.
        pub fn first(self: *Self) [*]FloatOf(T) {
            return @ptrCast(&self.buf[self.idx(0)]);
        }

        /// The lowest element, where an interface (and `rotm`) takes it.
        pub fn base(self: *Self) [*]T {
            return self.buf.ptr;
        }

        /// NaN in every element: the routine must not read them.
        pub fn poison(self: *Self) void {
            const nan = std.math.nan(FloatOf(T));
            for (0..self.n) |i| {
                self.was[i] = if (comptime isComplex(T)) .init(nan, nan) else nan;
                self.buf[self.idx(i)] = self.was[i];
            }
        }

        /// Every element as it was, bit for bit, and the gaps.
        pub fn unchanged(self: *const Self) !void {
            for (0..self.n) |i| if (!same(T, self.was[i], self.at(i), true)) {
                say("element {d} changed: ", .{i});
                return error.TestUnexpectedResult;
            };
            return self.gaps();
        }

        /// The gaps and the tail still hold their NaN.
        pub fn gaps(self: *const Self) !void {
            const s: usize = @abs(self.inc);
            for (self.buf, 0..) |v, k| {
                if (k % s == 0 and k / s < self.n) continue;
                if (!poisoned(T, v)) {
                    say("gap at {d} written: ", .{k});
                    return error.TestUnexpectedResult;
                }
            }
        }
    };
}

/// A case's patterns for an operand of `len`: `patterns` natively, four on
/// the rigs, where emulation is slow.
pub fn patternsFor(len: usize, buf: *[13]Pattern) []const Pattern {
    if (!rig()) return patterns(len, buf);
    buf[0..4].* = .{ .{ .one = .{ .kind = .nan, .at = len / 2 } }, .{ .one = .{ .kind = .inf, .at = len - 1 } }, .mixed, .denormal };
    return buf[0..4];
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
pub fn case(e: anyerror, routine: []const u8, comptime T: type, n: usize, inc: anytype) anyerror {
    say("{s} {s} n={d} inc={d},{d}: ", .{ routine, @typeName(T), n, inc[0], inc[1] });
    return e;
}

/// The core `eachCore` is running, for checks that know a core's quirks.
pub var core: ?cpu.CoreType = null;

/// Runs `check` on a vtable holding each linked core's kernels over the
/// portable ones, for every core this machine runs; natively, then on one
/// holding every `@Vector` kernel (`vector.zig`), which arm64 cores install,
/// on this machine's vectors. A failing core does not stop the rest, so one
/// run names them all.
pub fn eachCore(comptime check: fn (*const Vtable) anyerror!void) !void {
    var failed: ?anyerror = null;
    inline for (@typeInfo(cpu.CoreType).@"enum".field_names) |name| {
        const c = @field(cpu.CoreType, name);
        if (registry.isLinked(c) and (!@hasDecl(cpu.arch, "canRun") or cpu.arch.canRun(c))) {
            // on the rigs, where a core can take minutes, so a hang shows where
            if (rig()) say("[{s}] ", .{name});
            withCore(c, check) catch |e| {
                say("core {s}\n", .{name});
                failed = e;
            };
        }
    }
    if (!rig()) withVector(check) catch |e| {
        say("the @Vector kernels\n", .{});
        failed = e;
    };
    if (failed) |e| return e;
}

fn withVector(comptime check: fn (*const Vtable) anyerror!void) anyerror!void {
    var vt: Vtable = .{};
    portable.install(&vt);
    vector.installAll(&vt);
    portable.finish(&vt);
    return check(&vt);
}

/// Whether the tests run on a rig: the freestanding ones, and the Linux rigs
/// and remote Macs, which set `RSTD_TEST_PROGRESS`. `eachCore` names each
/// core as it starts there, and the biggest cases are left to native runs.
pub fn rig() bool {
    return switch (sysinfo.os.tag) {
        .freestanding => true,
        .linux, .macos => std.testing.environ.getPosix("RSTD_TEST_PROGRESS") != null,
        else => false,
    };
}

/// One core per call: inlined into `eachCore`'s loop, each core's table and
/// install would keep a frame slot of their own (~400 KB on i386 debug).
fn withCore(comptime c: cpu.CoreType, comptime check: fn (*const Vtable) anyerror!void) anyerror!void {
    core = c;
    defer core = null;
    var vt: Vtable = .{};
    portable.install(&vt);
    cpu.kernel.install(&vt, c);
    portable.finish(&vt);
    return check(&vt);
}

// ---- the extremes --------------------------------------------------------------

/// `kern` (a `value` or `index` kernel of the extremes, by its return type)
/// against a plain loop; `T` is the element type, real or complex.
pub fn checkExtreme(comptime T: type, comptime find: extremum.Find, routine: []const u8, kern: anytype) !void {
    const F = FloatOf(T);
    const complex = comptime isComplex(T);
    for (sizes) |n| {
        for ([_]usize{ 1, 2, 3 }) |inc| {
            var x: [capacity]T = undefined;
            fill(T, &x, n + 1);
            const W = Wide(F);
            var best: usize = 0;
            var best_key: W = undefined;
            for (0..n) |i| {
                const v = x[i * inc];
                const k: W = if (comptime complex) @abs(@as(W, v.re)) + @abs(@as(W, v.im)) else if (find.abs) @abs(@as(W, v)) else v;
                if (i == 0 or (if (find.max) k > best_key else k < best_key)) {
                    best = i;
                    best_key = k;
                }
            }
            const got = kern(@intCast(n), @as([*]F, @ptrCast(&x)), @intCast(inc));
            const ok = if (@TypeOf(got) == types.BlasLong)
                expectClose(f64, @floatFromInt(best + 1), @floatFromInt(got), 0)
            else
                expectCloseRef(F, best_key, got, tolerance(F, 1));
            ok catch |e| return case(e, routine, T, n, .{ inc, 0 });
        }
    }
}

test "the portable extremes at the extended type match a plain loop" {
    const E = Ext;
    inline for (.{ false, true }) |abs| inline for (.{ false, true }) |max| {
        const find: extremum.Find = .{ .abs = abs, .max = max };
        try checkExtreme(E, find, "value", &extremum.Extremum(E, 1, find).value);
        try checkExtreme(E, find, "index", &extremum.Extremum(E, 1, find).index);
        if (abs) {
            try checkExtreme(std.math.Complex(E), find, "value", &extremum.Extremum(E, 2, find).value);
            try checkExtreme(std.math.Complex(E), find, "index", &extremum.Extremum(E, 2, find).index);
        }
    };
}

test "tolerance goes by width, so a 64-bit c_longdouble gets f64's" {
    const L = c_longdouble;
    const want: L = if (@typeInfo(L).float.bits == 64) 1e-14 * 16 else 8 * std.math.floatEps(L) * 16;
    try zigstd.testExpectEqual(want, tolerance(L, 1));
    try zigstd.testExpectEqual(8 * std.math.floatEps(f80) * 16, tolerance(f80, 1));
}
