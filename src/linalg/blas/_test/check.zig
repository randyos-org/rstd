//! Test helpers shared by every level: type lists and shapes, random operands,
//! tolerant comparison, and running a check against every core's kernels.

const std = @import("std");
const sysinfo = @import("builtin");
const zigstd = @import("../../../zigstd.zig");
const cpu = @import("../cpu/__root__.zig");
const registry = @import("../registry.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const vector = @import("../vector.zig");
const extremum = @import("../blas1/extremum.zig");
const root = @import("root");
const bitwise = @import("bitwise.zig");
const pattern = @import("pattern.zig");
const reference = @import("reference.zig");

/// The `Io` the tests pass the routines: serial (`Io.failing`), which every
/// target has, rigs included.
pub const io = std.Io.failing;

/// The extended float the `q`/`x` portable kernels are tested at: `XDouble`
/// where it is one, else `f80`, whose kernels are generic Zig ones
/// (`types.kernelCc`) while the vtable's `q`/`x` slots are `f64`.
pub const Ext = if (types.isExtended(types.XDouble)) types.XDouble else f80;

/// `XDouble` and `ComplexXD` where they are types of their own, for the type
/// lists of tests that go through the vtable; empty where they alias `f64`.
pub const x_real = if (types.XDouble == f64) .{} else .{types.XDouble};
pub const x_complex = if (types.XDouble == f64) .{} else .{types.ComplexXD};

/// `f80` and its `Complex` where they have no slots, for the portable
/// kernels' type lists: Zig's calling convention, soft-float off x86.
pub const f80_unslotted = if (types.hasSlots(f80)) .{} else .{ f80, std.math.Complex(f80) };

/// The vtable's single and double types, real and complex.
pub const sd_types = .{ f32, f64, types.Complex32, types.Complex64 };
/// `sd_types` and the extended types with slots of their own.
pub const all_types = sd_types ++ x_real ++ x_complex;
/// `sd_types` by their precision letters.
/// The level-1 reductions the edge and denormal tests sweep.
pub const Reduce = enum {
    asum,
    sum,
    nrm2,

    pub fn op(r: Reduce) portable.Op {
        return switch (r) {
            .asum => .asum_k,
            .sum => .sum_k,
            .nrm2 => .nrm2_k,
        };
    }
};

/// The level-1 extremes: the value op, the index op, and what they find.
pub const extremes = [_]struct { portable.Op, portable.Op, extremum.Find }{
    .{ .amax_k, .iamax_k, .{ .abs = true, .max = true } },
    .{ .amin_k, .iamin_k, .{ .abs = true, .max = false } },
    .{ .max_k, .imax_k, .{ .abs = false, .max = true } },
    .{ .min_k, .imin_k, .{ .abs = false, .max = false } },
};

pub const sd_letters = .{ .{ "s", f32 }, .{ "d", f64 }, .{ "c", types.Complex32 }, .{ "z", types.Complex64 } };

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

/// Mixed magnitudes, from `seed`, so a sum's bits follow its order.
pub fn fillMixed(comptime T: type, x: []T, seed: u64) void {
    var prng: std.Random.DefaultPrng = .init(seed);
    const rand = prng.random();
    const F = FloatOf(T);
    for (x) |*v| {
        const re = std.math.ldexp(rand.float(F) * 2 - 1, rand.intRangeAtMost(i32, -12, 12));
        v.* = if (comptime isComplex(T)) .init(re, std.math.ldexp(rand.float(F) * 2 - 1, rand.intRangeAtMost(i32, -12, 12))) else re;
    }
}

/// Small integers in [-11, 11], from `seed`: products and short sums of them
/// are exact.
pub fn fillSeq(comptime F: type, buf: []F, seed: usize) void {
    for (buf, 0..) |*v, i| v.* = @floatFromInt(@as(i64, @intCast((i * 7 + seed * 13) % 23)) - 11);
}

/// `re + im i` as a `T`; `re` for a real `T`.
pub fn val(comptime T: type, re: f64, im: f64) T {
    return if (comptime isComplex(T)) .init(@floatCast(re), @floatCast(im)) else @floatCast(re);
}

fn random(comptime F: type, r: std.Random) F {
    if (comptime !types.isExtended(F)) return @floatCast(r.float(f64) * 4 - 2);
    // [1, 2) with every fraction bit random
    const bits = std.math.floatFractionalBits(F);
    const frac = r.int(u128) & ((@as(u128, 1) << bits) - 1);
    const v = 1 + @as(F, @floatFromInt(frac)) * comptime std.math.ldexp(@as(F, 1), -bits);
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

/// `expectClose` of `got` against a reference summed in `reference.Wide(F)`.
pub fn expectCloseRef(comptime F: type, want: reference.Wide(F), got: F, tol: reference.Wide(F)) !void {
    return close(@as(F, @floatCast(want)), got, @as(F, @floatCast(tol)));
}

pub fn expectCloseSlices(comptime T: type, want: []const T, got: []const T, tol: FloatOf(T)) !void {
    try zigstd.testExpectEqual(want.len, got.len);
    for (want, got, 0..) |w, g, i| expectClose(T, w, g, tol) catch |e| {
        say("element {d}: ", .{i});
        return e;
    };
}

/// Past a vector's last element, NaN that must stay as it was.
pub const tail = 4;

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
        pub fn init(al: std.mem.Allocator, n: usize, inc: isize, seed: u64, p: ?pattern.Pattern) !Self {
            const buf = try al.alloc(T, (n - 1) * @abs(inc) + 1 + tail);
            errdefer al.free(buf);
            var self: Self = .{ .buf = buf, .was = try al.alloc(T, n), .n = n, .inc = inc };
            @memset(std.mem.sliceAsBytes(self.buf), 0xff);
            fill(T, self.was, seed);
            if (p) |q| for (self.was, 0..) |*v, i| {
                v.* = pattern.shape(T, v.*, i, q);
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
        pub fn w(self: *const Self, i: usize) reference.Wide(T) {
            return reference.wide(T, self.was[i]);
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
            for (0..self.n) |i| if (!bitwise.same(T, self.was[i], self.at(i), true)) {
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
                if (!bitwise.poisoned(T, v)) {
                    say("gap at {d} written: ", .{k});
                    return error.TestUnexpectedResult;
                }
            }
        }
    };
}

/// The float types beyond the vtable's to test the interface at. `f80` on
/// every target: x87 on x86, soft-float elsewhere, and in slots only where it
/// is `XDouble`.
pub const any_floats = .{ f16, f80, f128 };

/// `util.extendedSum` written out plainly: accumulator by position, then
/// folded.
pub fn plainExtendedSum(comptime F: type, w: usize, term: util.Term, n: usize, x: []const F, ix: usize) F {
    var a: [4]F = @splat(0);
    const main = n / (8 / w) * (8 / w);
    for (0..n) |i| for (0..w) |c| {
        const v = x[i * ix * w + c];
        a[if (i < main) (i * w + c) % 4 else c] += if (term == .abs) @abs(v) else v * v;
    };
    return a[3] + ((a[2] + a[0]) + a[1]);
}

/// Says which case of a check failed: on stderr, or a qemu rig's console.
pub fn say(comptime f: []const u8, args: anytype) void {
    if (comptime sysinfo.target.os.tag != .freestanding) {
        std.debug.print(f, args);
    } else if (comptime @hasDecl(root, "testWrite")) {
        var buf: [256]u8 = undefined;
        root.testWrite(std.mem.print(&buf, f, args) catch return);
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
        if (registry.isLinked(c) and cpu.arch.canRun(c)) {
            // on the rigs, where a core can take minutes, so a hang shows where
            if (rig()) say("[{s}] ", .{name});
            withCore(c, check) catch |e| {
                say("core {s}: {t}\n", .{ name, e });
                failed = e;
            };
        }
    }
    if (!rig()) withVector(check) catch |e| {
        say("the @Vector kernels: {t}\n", .{e});
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

/// Whether the tests run on a rig: the freestanding and UEFI ones, WASI, and
/// the Linux rigs and remote Macs, which set `RSTD_TEST_PROGRESS`. `eachCore`
/// names each core as it starts there, and the biggest cases are left to
/// native runs.
pub fn rig() bool {
    return switch (sysinfo.target.os.tag) {
        .freestanding, .uefi, .wasi => true,
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
            const W = reference.Wide(F);
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
