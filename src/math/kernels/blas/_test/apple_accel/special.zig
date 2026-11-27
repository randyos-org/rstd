//! Accelerate against rstd's kernels on special values, strides and sizes, and the level-1 routing table.

const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../../__root__.zig");
const l3 = @import("../../blas3/core.zig");
const l2 = @import("../../blas2/core.zig");
const l1 = @import("../../blas1/interface.zig");
const accel = @import("../../apple_accel.zig");
const zigstd = @import("../../../../../zigstd.zig");
const check = @import("../check.zig");
const pattern = @import("../pattern.zig");
const reference = @import("../reference.zig");
const rt = @import("routing.zig");

const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const isComplex = types.isComplex;
const Float = types.Float;
const L1 = accel.L1;
const calls = &accel.calls;
const route = accel.route;
const routing = accel.routing;
const sent = accel.sent;
const x86_64 = accel.x86_64;
const Probe = rt.Probe;
const expectCalled = rt.expectCalled;
const val = rt.val;

// ---- Accelerate against rstd's kernels on special values ----------------------

/// The routed level-1 routines.
pub const R1 = enum { copy, swap, scal, dot, dotc, dsdot, nrm2, asum, iamax, rot, rotm, axpy, axpby };

/// `alpha`, `beta` (`scal`, `axpy`, `axpby`), `c`, `s` (`rot`) and `rotm`'s
/// flag, by variant.
pub const scalars1 = [_]struct { a: [2]f64, b: [2]f64, cs: [2]f64, flag: f64 }{
    .{ .a = .{ 1.5, -0.5 }, .b = .{ 0.5, 0.25 }, .cs = .{ 0.6, 0.8 }, .flag = -1 },
    .{ .a = .{ 0, 0 }, .b = .{ 0.5, 0 }, .cs = .{ 0, 1 }, .flag = 0 },
    .{ .a = .{ 1.5, 0 }, .b = .{ 0, 0 }, .cs = .{ 1, 0 }, .flag = 1 },
    .{ .a = .{ 0, 0 }, .b = .{ 0, 0 }, .cs = .{ -0.8, 0.6 }, .flag = -2 },
    .{ .a = .{ 1, 0 }, .b = .{ 1, 0 }, .cs = .{ 0.6, -0.8 }, .flag = -1 },
    .{ .a = .{ 0, 1 }, .b = .{ 0, -0.5 }, .cs = .{ 0.8, 0.6 }, .flag = 0 },
};

pub fn variants1(r: R1) usize {
    return switch (r) {
        .scal, .rot, .rotm, .axpy, .axpby => scalars1.len,
        else => 1,
    };
}

pub fn has1(comptime T: type, r: R1) bool {
    return switch (r) {
        .dotc => isComplex(T),
        .dsdot => T == f32,
        .rotm => !isComplex(T),
        else => true,
    };
}

/// A level-1 run's value or index.
pub const Out1 = struct { v: [2]f128 = .{ 0, 0 }, i: isize = 0 };

pub fn out1(v: anytype) Out1 {
    if (comptime isComplex(@TypeOf(v))) return .{ .v = .{ v.re, v.im } };
    return .{ .v = .{ v, 0 } };
}

pub fn run1(comptime T: type, comptime r: R1, k: usize, n: usize, x: *check.Vector(T), y: *check.Vector(T)) !Out1 {
    const io = check.io;
    const ni: isize = @intCast(n);
    const s = scalars1[k];
    const F = Float(T);
    const alpha = val(T, @floatCast(s.a[0]), @floatCast(s.a[1]));
    const beta = val(T, @floatCast(s.b[0]), @floatCast(s.b[1]));
    const xb = x.base();
    const yb = y.base();
    switch (r) {
        .copy => try l1.copy(io, T, ni, xb, x.inc, yb, y.inc),
        .swap => try l1.swap(io, T, ni, xb, x.inc, yb, y.inc),
        .scal => try l1.scal(io, T, ni, alpha, xb, x.inc),
        .dot => return out1(try l1.dot(io, T, ni, xb, x.inc, yb, y.inc)),
        .dotc => if (comptime isComplex(T)) return out1(try l1.dotc(io, T, ni, xb, x.inc, yb, y.inc)),
        .dsdot => if (comptime T == f32) return out1(try l1.dsdot(io, ni, xb, x.inc, yb, y.inc)),
        .nrm2 => return out1(try l1.nrm2(io, T, ni, xb, x.inc)),
        .asum => return out1(try l1.asum(io, T, ni, xb, x.inc)),
        .iamax => return .{ .i = try l1.index(io, T, .largest_magnitude, ni, xb, x.inc) },
        .rot => try l1.rot(io, T, ni, xb, x.inc, yb, y.inc, @floatCast(s.cs[0]), @floatCast(s.cs[1])),
        .rotm => if (comptime !isComplex(T)) {
            const param = [5]F{ @floatCast(s.flag), 0.75, -0.5, 0.25, 1.25 };
            try l1.rotm(io, T, ni, xb, x.inc, yb, y.inc, &param);
        },
        .axpy => try l1.axpy(io, T, ni, alpha, xb, x.inc, yb, y.inc),
        .axpby => try l1.axpby(io, T, ni, alpha, xb, x.inc, beta, yb, y.inc),
    }
    return .{};
}

/// `a` and `b` alike: NaN for NaN, the same infinity, else within `ulps` of
/// `scale`; with `lax`, any non-finite for any (complex infinities meet in
/// orders of their own).
pub fn alike(comptime F: type, a: f128, b: f128, ulps: f128, scale: f128, lax: bool) bool {
    if (!std.math.isFinite(a) or !std.math.isFinite(b)) {
        if (lax and !std.math.isFinite(a) and !std.math.isFinite(b)) return true;
        if (std.math.isNan(a) or std.math.isNan(b)) return std.math.isNan(a) and std.math.isNan(b);
        return a == b;
    }
    return @abs(a - b) <= ulps * (std.math.floatEps(F) * scale + std.math.floatTrueMin(F));
}

/// Element `i` of `v` as the routine read it (a zero stride reads one).
pub fn input(comptime T: type, v: *const check.Vector(T), i: usize) T {
    return if (v.inc == 0) v.was[v.n - 1] else v.was[i];
}

/// The case's results against rstd's: values within their accuracy, an
/// index exactly, vectors element by element, gaps untouched.
pub fn compare1(comptime T: type, comptime r: R1, n: usize, want: Out1, got: Out1, x1: *const check.Vector(T), y1: *const check.Vector(T), x2: *const check.Vector(T), y2: *const check.Vector(T)) !void {
    const F = Float(T);
    const cx = comptime isComplex(T);
    const nf: f128 = @floatFromInt(n);
    switch (r) {
        .iamax => if (want.i != got.i) {
            check.say("index {d}, rstd {d}: ", .{ got.i, want.i });
            return error.TestUnexpectedResult;
        },
        .dot, .dotc, .dsdot, .nrm2, .asum => {
            var scale: f128 = 0;
            for (0..n) |i| scale += switch (r) {
                .nrm2, .asum => reference.mag(input(T, x1, i)),
                else => reference.mag(input(T, x1, i)) * reference.mag(input(T, y1, i)),
            };
            if (r == .nrm2) scale = reference.mag(want.v[0]);
            // `sdot`, `snrm2` and `scnrm2` sum in `f64` and round once, as rstd does
            const single = (T == f32 and r == .dot) or (F == f32 and r == .nrm2);
            const ulps: f128 = if (single) 2 else nf + 4;
            const scale2 = if (single) reference.mag(want.v[0]) + scale * 0x1p-40 else scale;
            const RF = if (r == .dsdot) f64 else F;
            if (!alike(RF, want.v[0], got.v[0], ulps, scale2, cx) or !alike(RF, want.v[1], got.v[1], ulps, scale2, cx)) {
                check.say("{e},{e}, rstd {e},{e}: ", .{ got.v[0], got.v[1], want.v[0], want.v[1] });
                return error.TestUnexpectedResult;
            }
        },
        else => {},
    }
    for ([_][2]*const check.Vector(T){ .{ x1, x2 }, .{ y1, y2 } }, 0..) |pair, which| {
        if (pair[0].inc == 0) continue;
        for (0..n) |i| {
            const scale = 4 * (reference.mag(input(T, x1, i)) + reference.mag(input(T, y1, i)));
            const w = pair[0].at(i);
            const g = pair[1].at(i);
            const ok = if (comptime cx)
                alike(F, w.re, g.re, 8, scale, true) and alike(F, w.im, g.im, 8, scale, true)
            else
                alike(F, w, g, 8, scale, false);
            if (!ok) {
                check.say("{s}[{d}] {any}, rstd {any}: ", .{ if (which == 0) "x" else "y", i, g, w });
                return error.TestUnexpectedResult;
            }
        }
        try pair[1].gaps();
    }
}

/// Room for four vectors of 1025 at stride 3, and their copies as they were.
pub var mem1: [4 * (2 * 1025 * 3 + 64) * @sizeOf(Complex64)]u8 align(16) = undefined;

/// One case through rstd's kernels and through Accelerate, compared.
pub fn case1(comptime T: type, comptime r: R1, k: usize, n: usize, inc: [2]isize, p: ?pattern.Pattern, on_y: bool) !void {
    var fba: std.heap.FixedBufferAllocator = .init(&mem1);
    const al = fba.allocator();
    const px = if (on_y) null else p;
    const py = if (on_y) p else null;
    var x1: check.Vector(T) = try .init(al, n, inc[0], 7 * n + 1, px);
    var y1: check.Vector(T) = try .init(al, n, inc[1], 7 * n + 2, py);
    var x2: check.Vector(T) = try .init(al, n, inc[0], 7 * n + 1, px);
    var y2: check.Vector(T) = try .init(al, n, inc[1], 7 * n + 2, py);
    _ = route(false);
    const want = try run1(T, r, k, n, &x1, &y1);
    _ = route(true);
    const got = try run1(T, r, k, n, &x2, &y2);
    try compare1(T, r, n, want, got, &x1, &y1, &x2, &y2);
}

/// Distinct failures seen, so one run reports each once.
pub const Seen = struct {
    keys: [256][96]u8 = undefined,
    lens: [256]usize = undefined,
    len: usize = 0,
    failed: bool = false,

    /// Records `key`; true the first time.
    pub fn first(s: *Seen, key: []const u8) bool {
        s.failed = true;
        for (s.keys[0..s.len], s.lens[0..s.len]) |*k, l| if (std.mem.eql(u8, k[0..l], key)) return false;
        if (s.len == s.keys.len) return false;
        const l = @min(key.len, 96);
        @memcpy(s.keys[s.len][0..l], key[0..l]);
        s.lens[s.len] = l;
        s.len += 1;
        return true;
    }
};

pub fn patternName(p: ?pattern.Pattern) []const u8 {
    const q = p orelse return "plain";
    return switch (q) {
        .one => |o| @tagName(o.kind),
        else => @tagName(q),
    };
}

/// Every size, stride and pattern of `r` over `T`; a failure is reported once
/// per variant, pattern kind and stride sign.
pub fn sweep1(comptime T: type, comptime r: R1, seen: *Seen) void {
    const two = switch (r) {
        .nrm2, .asum, .iamax, .scal => false,
        else => true,
    };
    const strides = [_][2]isize{ .{ 1, 1 }, .{ 2, 3 }, .{ -1, -1 }, .{ -3, 2 }, .{ 2, -3 }, .{ 0, 1 }, .{ 1, 0 } };
    const reads_x = r == .dot or r == .dotc or r == .dsdot or r == .copy or r == .axpy or r == .axpby;
    const reads_y = r == .dot or r == .dotc or r == .dsdot;
    for (0..variants1(r)) |k| for (strides) |inc| {
        // zero strides only where the vector is read, not written
        if (inc[0] == 0 and !reads_x) continue;
        if (inc[1] == 0 and !reads_y) continue;
        if (!two and inc[1] != 1) continue;
        for ([_]usize{ 1, 2, 5, 16, 33, 67, 300, 1025 }) |n| {
            var buf: [13]pattern.Pattern = undefined;
            const ps = pattern.patterns(n, &buf);
            for (0..ps.len + 1) |pi| for ([_]bool{ false, true }) |on_y| {
                const p: ?pattern.Pattern = if (pi == ps.len) null else ps[pi];
                if (on_y and (!two or p == null)) continue;
                case1(T, r, k, n, inc, p, on_y) catch {
                    var kb: [96]u8 = undefined;
                    const key = std.mem.print(&kb, "{s} {s} v{d} {s} in {s} inc {d},{d}", .{ @tagName(r), @typeName(T), k, patternName(p), if (on_y) "y" else "x", std.math.sign(inc[0]), std.math.sign(inc[1]) }) catch "?";
                    if (seen.first(key)) check.say("n={d} inc={d},{d}: {s}\n", .{ n, inc[0], inc[1], key });
                };
            };
        }
    };
}

test "Accelerate's level-1 results match rstd's kernels on special values, strides and sizes" {
    _ = blas.vtable();
    const was = routing();
    defer _ = route(was);
    if (!route(true)) return error.SkipZigTest;
    var seen: Seen = .{};
    inline for (.{ f32, f64, Complex32, Complex64 }) |T| {
        inline for (comptime std.enums.values(R1)) |r| if (comptime has1(T, r)) sweep1(T, r, &seen);
    }
    if (seen.failed) return error.TestUnexpectedResult;
}

/// Whether a level-1 call at stride `inc` stays off Accelerate whatever
/// `sends` says: `cblas_dsdot` under Rosetta 2, the strides `nrm2`, `asum`
/// and `scal` (a quick return) leave out, and `.plain`'s single-precision
/// `asum` and complex `dot`, summed in rstd's `f64` slots (issue #584).
pub fn kept1(comptime T: type, r: R1, inc: isize) bool {
    return switch (r) {
        .dot => (T == f32 and accel.rosetta()) or T == Complex32,
        .dsdot => T == f32 and accel.rosetta(),
        .dotc => T == Complex32,
        .nrm2 => inc <= 0 or (Float(T) == f32 and accel.rosetta()),
        .asum => inc <= 0 or Float(T) == f32,
        .scal => inc <= 0,
        .iamax => true,
        else => false,
    };
}

/// Runs `r` over `T` on 6 elements at stride `inc` in both vectors; fails
/// unless it called Accelerate exactly when `want`.
pub fn expect1(comptime T: type, comptime r: R1, inc: isize, want: bool) !void {
    var fba: std.heap.FixedBufferAllocator = .init(&mem1);
    const al = fba.allocator();
    var x: check.Vector(T) = try .init(al, 6, inc, 1, null);
    var y: check.Vector(T) = try .init(al, 6, inc, 2, null);
    const before = calls.load(.monotonic);
    _ = try run1(T, r, 0, 6, &x, &y);
    expectCalled(want, before, @tagName(r)) catch |e| {
        check.say("{s} at stride {d}\n", .{ @typeName(T), inc });
        return e;
    };
}

test "level-1 calls go to Accelerate by routine, letter and stride as `sends` has them" {
    _ = blas.vtable();
    const was = routing();
    defer _ = route(was);
    if (!route(true)) return error.SkipZigTest;
    inline for (.{ f32, f64, Complex32, Complex64 }) |T| {
        inline for (comptime std.enums.values(L1)) |l| {
            const r = @field(R1, @tagName(l));
            if (comptime has1(T, r)) for ([_]isize{ 1, 2, -1 }) |inc| {
                try expect1(T, r, inc, sent(l, T, inc == 1) and !kept1(T, r, inc));
            };
        }
    }
    if (x86_64) return;
    // the arm64 flyoff's rows behind `sends`, either way
    try expect1(Complex64, .nrm2, 1, false);
    try expect1(f64, .nrm2, 1, false);
    try expect1(f32, .dot, 1, false);
    try expect1(f32, .dsdot, 1, true);
    try expect1(f64, .dot, 2, true);
    try expect1(Complex64, .dot, 1, true);
    try expect1(Complex64, .dot, 2, false);
    try expect1(f64, .copy, 2, false);
    try expect1(Complex32, .copy, 2, true);
    try expect1(f32, .copy, 1, true);
    try expect1(f32, .swap, 2, false);
    try expect1(Complex64, .swap, 2, false);
    try expect1(f64, .axpy, 2, true);
    try expect1(Complex32, .axpy, 2, false);
    try expect1(Complex32, .rot, 2, false);
    try expect1(Complex64, .rot, 2, true);
    try expect1(f64, .scal, 2, true);
    // `.plain`'s, in rstd's `f64` slot (issue #584)
    try expect1(f32, .asum, 2, false);
}

test "quirk cases run rstd's kernels, and the Hermitian diagonal comes out real" {
    _ = blas.vtable();
    if (!routing()) return error.SkipZigTest;
    const al = zigstd.testAllocator;
    const io = check.io;
    const inf = std.math.inf(f64);
    var z: [4]Complex64 = .{ .init(inf, 1), .init(1, 2), .init(-1, 0.5), .init(2, -1) };
    var w: [4]Complex64 = .{ .init(1, 1), .init(1, -1), .init(0.5, 0.5), .init(-2, 1) };
    var x = [_]f64{ 1, inf, -1, 2 };
    var y = [_]f64{ 1, 0, 2, 0.5 };
    var a: [16]f64 = @splat(1);
    for (0..4) |i| a[i * 5] = 2;
    var b: [16]f64 = @splat(0.5);
    b[5] = inf;
    var p: Probe = .init(true);

    // the zero-part quirks and real `trmv`'s show on arm64 only (issue #388)
    try l1.scal(io, Complex64, 4, .init(1.5, 0), &z, 1);
    try p.saw(x86_64 and sent(.scal, Complex64, true), "complex scal, alpha with a zero part");
    try zigstd.testExpectTrue(std.math.isNan(z[0].im));
    try l1.axpy(io, Complex64, 4, .init(0, 1), &z, 1, &w, 1);
    try p.saw(x86_64 and sent(.axpy, Complex64, true), "complex axpy, alpha with a zero part");
    try l1.axpby(io, f64, 4, 0, &x, 1, 2, &y, 1);
    try p.saw(false, "axpby, zero alpha");
    try zigstd.testExpectEqualSlices(f64, &.{ 2, 0, 4, 1 }, &y);
    try l1.axpby(io, Complex64, 4, .init(1, 1), &z, 1, .init(2, 0), &w, 1);
    try p.saw(x86_64, "complex axpby, beta with a zero part");
    // real `symv`'s too (issue #410): arm64 gives NaN where this is inf
    var sv: [4]f64 = @splat(0);
    try l2.symv(al, io, f64, .u, 4, 1, &a, 4, &x, 1, 0, &sv, 1);
    try p.saw(x86_64 and accel.routes(.symv, f64), "symv, x not finite");
    for (sv) |v| try zigstd.testExpectEqual(inf, v);
    try l2.trmv(al, io, f64, .u, .n, .n, 4, &a, 4, &x, 1);
    try p.saw(x86_64, "trmv, x not finite");
    try l2.ger(al, io, f64, 4, 4, 1, &x, 1, &y, 1, &a, 4);
    try p.saw(false, "ger, x not finite");
    // a zero factor beside an infinity, through rstd; without one, Accelerate
    var sx = [_]f64{ 0, inf, -1, 2 };
    @memset(&a, 1);
    try l2.syr(al, io, f64, .l, 4, 1, &sx, 1, &a, 4);
    try p.saw(false, "syr, x(0) zero and x not finite");
    try zigstd.testExpectEqual(@as(f64, 1), a[1]);
    // `x(3)` and `y(3)` are zero, `x` at a negative stride
    var sy = [_]f64{ 1, 2, 3, -0.0 };
    try l2.syr2(al, io, f64, .u, 4, 1, &sx, -1, &sy, 1, &a, 4);
    try p.saw(false, "syr2, x(j) and y(j) zero and x not finite");
    sy[3] = 1;
    try l2.syr2(al, io, f64, .u, 4, 1, &sx, -1, &sy, 1, &a, 4);
    try p.saw(true, "syr2, no j with both zero");
    try l3.trmm(al, io, f64, .l, .l, .n, .n, 4, 4, 1, &a, 4, &b, 4);
    try p.saw(false, "trmm, B not finite");
    // the same cases, finite, call it
    x = .{ 1, 3, -1, 2 };
    @memset(&a, 1);
    try l2.trmv(al, io, f64, .u, .n, .n, 4, &a, 4, &x, 1);
    try p.saw(true, "trmv");
    @memset(&b, 0.5);
    try l3.trmm(al, io, f64, .l, .l, .n, .n, 4, 4, 1, &a, 4, &b, 4);
    try p.saw(true, "trmm");

    // `herk` and `her2k` make the diagonal real; the reference never reads it
    var h: [12]Complex64 = undefined;
    var c: [9]Complex64 = undefined;
    for ([_]bool{ false, true }) |two| for ([_]bool{ false, true }) |special| {
        check.fill(Complex64, &h, 3);
        if (special) h[4] = .init(inf, 1);
        check.fill(Complex64, &c, 4);
        for (0..3) |i| c[i * 4].im = 3;
        if (two) {
            try l3.her2k(al, io, Complex64, .u, .n, 3, 4, .init(0.5, 0.25), &h, 3, &h, 3, 0.5, &c, 3);
        } else {
            try l3.herk(al, io, Complex64, .l, .n, 3, 4, 1.5, &h, 3, 0.5, &c, 3);
        }
        try p.saw(true, "herk and her2k");
        for (0..3) |i| try zigstd.testExpectEqual(@as(f64, 0), c[i * 4].im);
    };
    // and `her2`, packed too
    var hx: [3]Complex64 = .{ .init(inf, 1), .init(1, 2), .init(-1, 0.5) };
    check.fill(Complex64, &c, 5);
    try l2.her2(al, io, Complex64, .l, 3, .init(0.5, 0.25), &hx, 1, &hx, 1, &c, 3);
    try p.saw(true, "her2");
    for (0..3) |i| try zigstd.testExpectEqual(@as(f64, 0), c[i * 4].im);
    check.fill(Complex64, &c, 5);
    try l2.hpr2(al, io, Complex64, .u, 3, .init(0.5, 0.25), &hx, 1, &hx, 1, &c);
    try p.saw(true, "hpr2");
    for ([_]usize{ 0, 2, 5 }) |i| try zigstd.testExpectEqual(@as(f64, 0), c[i].im);
    // `beta` 1 and nothing added leaves `C` as it was, through rstd
    check.fill(Complex64, &c, 4);
    for (0..3) |i| c[i * 4].im = 3;
    try l3.herk(al, io, Complex64, .l, .n, 3, 0, 1.5, &h, 3, 1, &c, 3);
    try p.saw(false, "herk, beta 1 and k 0");
    for (0..3) |i| try zigstd.testExpectEqual(@as(f64, 3), c[i * 4].im);
    // a non-finite `C` with a nonzero `beta`, through rstd; with a zero one, Accelerate
    check.fill(Complex64, &h, 3);
    check.fill(Complex64, &c, 4);
    c[1] = .init(1, inf);
    try l3.herk(al, io, Complex64, .l, .n, 3, 4, 1.5, &h, 3, 0.5, &c, 3);
    try p.saw(false, "herk, C not finite");
    try zigstd.testExpectTrue(!std.math.isNan(c[1].re));
    c[1] = .init(1, inf);
    try l3.herk(al, io, Complex64, .l, .n, 3, 4, 1.5, &h, 3, 0, &c, 3);
    try p.saw(true, "herk, beta 0");
}
