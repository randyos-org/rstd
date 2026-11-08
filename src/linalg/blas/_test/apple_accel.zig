//! `apple_accel.zig`'s tests: binding, routing, the split over rstd's
//! workers, and special values against rstd's kernels.

const std = @import("std");
const sysinfo = @import("builtin");
const dll = @import("../../../dll.zig");
const types = @import("../types.zig");
const blas = @import("../__root__.zig");
const threading = @import("../threading.zig");
const test_io = @import("io.zig");
const accel = @import("../apple_accel.zig");
const zigstd = @import("../../../zigstd.zig");
const check = @import("check.zig");
const pattern = @import("pattern.zig");
const reference = @import("reference.zig");

const Io = std.Io;

const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const isComplex = types.isComplex;
const Float = types.Float;

const AppleAccelerate = accel.AppleAccelerate;
const Blas = accel.Blas;
const L1 = accel.L1;
const Table = accel.Table;
const ThreadApi = accel.ThreadApi;
const Threading = accel.Threading;
const available = accel.available;
const bindBlas = accel.bindBlas;
const call = accel.call;
const calls = &accel.calls;
const fortran = accel.fortran;
const ilp64Names = accel.ilp64Names;
const load = accel.load;
const lp64Names = accel.lp64Names;
const matrix_quantum = accel.matrix_quantum;
const route = accel.route;
const routing = accel.routing;
const sent = accel.sent;
const span = accel.span;
const unpinned = &accel.unpinned;
const vector_quantum = accel.vector_quantum;
const x86_64 = accel.x86_64;

/// The integer type of `tag`'s table.
fn Int(comptime tag: std.meta.Tag(Blas)) type {
    return if (tag == .ilp64) i64 else c_int;
}

fn fill(comptime F: type, rng: std.Random, out: []F) void {
    for (out) |*v| v.* = rng.float(F) * 2 - 1;
}

test "load binds the CBLAS table on macOS, and nothing elsewhere" {
    const acc = load(.{}) orelse return zigstd.testExpectTrue(!available);
    switch (acc.blas) {
        inline else => |b| {
            try zigstd.testExpectTrue(b.fns.cblas_dgemm != null and b.fns.cblas_ddot != null);
            try zigstd.testExpectTrue(b.fns.cblas_zdotu_sub != null and b.fns.cblas_isamax != null);
        },
    }
    if (acc.threads.BLASSetThreading != null) try zigstd.testExpectEqual(Threading.single_threaded, acc.threading().?);
}

/// Accelerate's exports of `Table`'s names on one macOS, from the SDK's
/// `vecLib.framework` `libBLAS.tbd` and `libLAPACK.tbd` (MacOSX11.3, 12.3 and
/// 13.3.sdk, <https://github.com/alexey-lysiuk/macos-sdk> at 896cd40d; arm64
/// and x86_64 list the same). Every field but `cblas_i?amin` is there, each
/// `?potrf_` also as `?potrf`; 13.3 adds each with `$NEWLAPACK` and
/// `$NEWLAPACK$ILP64`, LAPACK's without the underscore. None has
/// `BLASSetThreading` (macOS 15+). Re-derive from a newer SDK's `.tbd`s
/// rather than editing by hand.
const RecordedAccelerate = struct {
    /// macOS 13.3 or later
    new: bool,

    fn has(r: RecordedAccelerate, name: []const u8) bool {
        const cut = std.mem.findScalar(u8, name, '$') orelse name.len;
        const base = name[0..cut];
        const suffix = name[cut..];
        const known = std.mem.eql(u8, suffix, "") or
            (r.new and (std.mem.eql(u8, suffix, "$NEWLAPACK") or std.mem.eql(u8, suffix, "$NEWLAPACK$ILP64")));
        if (!known or std.mem.endsWith(u8, base, "amin")) return false;
        inline for (@typeInfo(Table(i64)).@"struct".field_names) |f| {
            if (comptime fortran(f)) |b| {
                if (std.mem.eql(u8, base, b)) return true;
                if (std.mem.eql(u8, base, f)) return suffix.len == 0;
            } else if (std.mem.eql(u8, base, f)) return true;
        }
        return false;
    }

    /// `Dll.lookup`: a pointer no test calls, for a name this macOS has.
    pub fn lookup(r: *const RecordedAccelerate, comptime T: type, name: [:0]const u8) ?T {
        return if (r.has(name)) @ptrFromInt(0x1000) else null;
    }
};

test "macOS 11 and 12 bind the classic 32-bit names" {
    const old: RecordedAccelerate = .{ .new = false };
    for ([_]bool{ false, true }) |lp64| {
        const b = bindBlas(&old, lp64) orelse return error.TestUnexpectedResult;
        try zigstd.testExpectTrue(b == .lp64);
        // 158 fields, less the four `i?amin`
        try zigstd.testExpectEqual(@as(usize, 154), b.lp64.count());
        try zigstd.testExpectEqualStrings("cblas_dgemm", b.lp64.name(lp64Names, .cblas_dgemm).?);
        try zigstd.testExpectEqualStrings("catlas_daxpby", b.lp64.name(lp64Names, .catlas_daxpby).?);
        try zigstd.testExpectEqualStrings("zpotrf_", b.lp64.name(lp64Names, .zpotrf_).?);
        try zigstd.testExpectTrue(b.lp64.fns.cblas_isamin == null);
    }
    // nor the thread API: `VECLIB_MAXIMUM_THREADS` alone keeps it on one thread
    const acc: AppleAccelerate = .{ .blas = bindBlas(&old, false).?, .threads = dll.bind(ThreadApi, &old, dll.exactName).fns };
    try zigstd.testExpectTrue(acc.threading() == null);
    try zigstd.testExpectError(error.Unsupported, acc.setThreading(.single_threaded));
}

test "macOS 13.3 binds ILP64, or the new 32-bit names" {
    const new: RecordedAccelerate = .{ .new = true };
    const ilp = bindBlas(&new, false) orelse return error.TestUnexpectedResult;
    try zigstd.testExpectTrue(ilp == .ilp64);
    try zigstd.testExpectEqual(@as(usize, 154), ilp.ilp64.count());
    try zigstd.testExpectEqualStrings("cblas_dgemm$NEWLAPACK$ILP64", ilp.ilp64.name(ilp64Names, .cblas_dgemm).?);
    try zigstd.testExpectEqualStrings("catlas_daxpby$NEWLAPACK$ILP64", ilp.ilp64.name(ilp64Names, .catlas_daxpby).?);
    try zigstd.testExpectEqualStrings("dpotrf$NEWLAPACK$ILP64", ilp.ilp64.name(ilp64Names, .dpotrf_).?);
    const lp = bindBlas(&new, true) orelse return error.TestUnexpectedResult;
    try zigstd.testExpectTrue(lp == .lp64);
    try zigstd.testExpectEqual(@as(usize, 154), lp.lp64.count());
    try zigstd.testExpectEqualStrings("cblas_dgemm$NEWLAPACK", lp.lp64.name(lp64Names, .cblas_dgemm).?);
    try zigstd.testExpectEqualStrings("spotrf$NEWLAPACK", lp.lp64.name(lp64Names, .spotrf_).?);
    // a library with none of it
    const Nothing = struct {
        pub fn lookup(_: *const @This(), comptime T: type, _: [:0]const u8) ?T {
            return null;
        }
    };
    try zigstd.testExpectTrue(bindBlas(&Nothing{}, false) == null);
}

test "the 32-bit table binds the new names first, and agrees with ILP64" {
    if (comptime !available) return error.SkipZigTest;
    const acc = load(.{}) orelse return error.SkipZigTest;
    if (acc.blas != .ilp64) return error.SkipZigTest;
    var lib = try dll.openDarwinFramework("Accelerate");
    const lp = dll.bind(Table(c_int), &lib, lp64Names);
    try zigstd.testExpectEqualStrings("cblas_ddot$NEWLAPACK", lp.name(lp64Names, .cblas_ddot).?);
    var x: [100]f64 = undefined;
    for (&x, 0..) |*v, i| v.* = @floatFromInt(i % 7);
    try zigstd.testExpectEqual(acc.blas.ilp64.fns.cblas_ddot.?(100, &x, 1, &x, 1), lp.fns.cblas_ddot.?(100, &x, 1, &x, 1));
}

test "ddot, dgemm, zdotu_sub and isamax agree with rstd" {
    const acc = load(.{}) orelse return error.SkipZigTest;
    // rstd's own kernels
    _ = blas.vtable();
    const was = routing();
    _ = route(false);
    defer _ = route(was);
    const gpa = zigstd.testAllocator;
    var prng: std.Random.DefaultPrng = .init(0x5eed);
    const rng = prng.random();
    switch (acc.blas) {
        inline else => |b, tag| {
            const I = Int(tag);
            const t = b.fns;

            // ddot, strided both ways
            var x: [600]f64 = undefined;
            var y: [600]f64 = undefined;
            fill(f64, rng, &x);
            fill(f64, rng, &y);
            const got = t.cblas_ddot.?(300, &x, 2, &y, -1);
            const want = try blas.dot(check.io, f64, 300, &x, 2, &y, -1);
            try zigstd.testExpectApproxEqRel(want, got, 1e-12);

            // zdotu_sub and zdotc_sub
            var zx: [257]Complex64 = undefined;
            var zy: [257]Complex64 = undefined;
            fill(f64, rng, @as(*[514]f64, @ptrCast(&zx)));
            fill(f64, rng, @as(*[514]f64, @ptrCast(&zy)));
            var zd: Complex64 = undefined;
            t.cblas_zdotu_sub.?(zx.len, &zx, 1, &zy, 1, &zd);
            const zw = try blas.dot(check.io, Complex64, zx.len, &zx, 1, &zy, 1);
            try zigstd.testExpectApproxEqAbs(zw.re, zd.re, 1e-12);
            try zigstd.testExpectApproxEqAbs(zw.im, zd.im, 1e-12);
            t.cblas_zdotc_sub.?(zx.len, &zx, 1, &zy, 1, &zd);
            const zc = try blas.dotc(check.io, Complex64, zx.len, &zx, 1, &zy, 1);
            try zigstd.testExpectApproxEqAbs(zc.re, zd.re, 1e-12);
            try zigstd.testExpectApproxEqAbs(zc.im, zd.im, 1e-12);

            // isamax: 0-based where rstd's is 1-based
            var sx: [1000]f32 = undefined;
            fill(f32, rng, &sx);
            sx[613] = -3;
            const ia: I = t.cblas_isamax.?(sx.len, &sx, 1);
            try zigstd.testExpectEqual(@as(I, 613), ia);
            try zigstd.testExpectEqual(try blas.index(check.io, f32, .amax, sx.len, &sx, 1), ia + 1);

            // dgemm
            const m = 37;
            const n = 29;
            const k = 41;
            const a = try gpa.alloc(f64, m * k);
            defer gpa.free(a);
            const bm = try gpa.alloc(f64, k * n);
            defer gpa.free(bm);
            const c0 = try gpa.alloc(f64, m * n);
            defer gpa.free(c0);
            const c1 = try gpa.alloc(f64, m * n);
            defer gpa.free(c1);
            fill(f64, rng, a);
            fill(f64, rng, bm);
            fill(f64, rng, c0);
            @memcpy(c1, c0);
            // A^T, stored k by m
            t.cblas_dgemm.?(.col_major, .trans, .no_trans, m, n, k, 1.5, a.ptr, k, bm.ptr, k, -0.5, c0.ptr, m);
            try blas.gemm(gpa, check.io, f64, .t, .n, m, n, k, 1.5, a.ptr, k, bm.ptr, k, -0.5, c1.ptr, m);
            for (c0, c1) |g, w| try zigstd.testExpectApproxEqAbs(w, g, 1e-12 * k);
        },
    }
}

/// Each routine's result bits.
const Digest = struct { ddot: u64, dnrm2: u64, dasum: u64, dgemv_n: u64, dgemv_t: u64, dgemm: u64 };

const ReproData = struct {
    x: []f64,
    y: []f64,
    a: []f64,
    b: []f64,
    c: []f64,
    v: []f64,
    const n1 = 1_000_000;
    const n2 = 2000;
};

fn hash(s: []const f64) u64 {
    return std.hash.Wyhash.hash(0, std.mem.sliceAsBytes(s));
}

fn digest(comptime I: type, t: *const Table(I), d: *const ReproData) Digest {
    const n1 = ReproData.n1;
    const n2 = ReproData.n2;
    var out: Digest = undefined;
    out.ddot = @bitCast(t.cblas_ddot.?(n1, d.x.ptr, 1, d.y.ptr, 1));
    out.dnrm2 = @bitCast(t.cblas_dnrm2.?(n1, d.x.ptr, 1));
    out.dasum = @bitCast(t.cblas_dasum.?(n1, d.x.ptr, 1));
    t.cblas_dgemv.?(.col_major, .no_trans, n2, n2, 1, d.a.ptr, n2, d.x.ptr, 1, 0, d.v.ptr, 1);
    out.dgemv_n = hash(d.v);
    t.cblas_dgemv.?(.col_major, .trans, n2, n2, 1, d.a.ptr, n2, d.x.ptr, 1, 0, d.v.ptr, 1);
    out.dgemv_t = hash(d.v);
    t.cblas_dgemm.?(.col_major, .no_trans, .no_trans, n2, n2, n2, 1, d.a.ptr, n2, d.b.ptr, n2, 0, d.c.ptr, n2);
    out.dgemm = hash(d.c);
    return out;
}

fn digestOnThread(comptime I: type, acc: *const AppleAccelerate, t: *const Table(I), d: *const ReproData, out: *Digest) void {
    acc.singleThreadHere();
    out.* = digest(I, t, d);
}

// Threaded, Accelerate's `gemv` with `A` untransposed rounds differently from
// one thread (macOS 27, M6); `ddot`, `dnrm2`, `dasum`, `A^T` `gemv` and `gemm`
// did not differ for any thread count tried. `load` pins one thread.
test "single-threaded results are bit-identical across runs and threads" {
    if (comptime !available or sysinfo.single_threaded) return error.SkipZigTest;
    const acc = load(.{}) orelse return error.SkipZigTest;
    const gpa = zigstd.testAllocator;
    var prng: std.Random.DefaultPrng = .init(0xacce1);
    const rng = prng.random();
    const n1 = ReproData.n1;
    const n2 = ReproData.n2;
    var d: ReproData = undefined;
    d.x = try gpa.alloc(f64, n1);
    defer gpa.free(d.x);
    d.y = try gpa.alloc(f64, n1);
    defer gpa.free(d.y);
    d.a = try gpa.alloc(f64, n2 * n2);
    defer gpa.free(d.a);
    d.b = try gpa.alloc(f64, n2 * n2);
    defer gpa.free(d.b);
    d.c = try gpa.alloc(f64, n2 * n2);
    defer gpa.free(d.c);
    d.v = try gpa.alloc(f64, n2);
    defer gpa.free(d.v);
    for ([_][]f64{ d.x, d.y, d.a, d.b }) |s| fill(f64, rng, s);

    switch (acc.blas) {
        inline else => |b, tag| {
            const I = Int(tag);
            const t = &b.fns;
            const first = digest(I, t, &d);
            try zigstd.testExpectEqual(first, digest(I, t, &d));
            var other: Digest = undefined;
            const th = try std.Thread.spawn(.{}, digestOnThread, .{ I, &acc, t, &d, &other });
            th.join();
            try zigstd.testExpectEqual(first, other);
            // the process-wide cap holds even where a thread asks for more
            if (acc.setThreading(.multi_threaded)) {
                defer acc.singleThreadHere();
                try zigstd.testExpectEqual(first, digest(I, t, &d));
            } else |_| {}
        },
    }
}

// ---- routing tests -------------------------------------------------------------

/// Fails unless Accelerate was called since `before` exactly when `want`.
fn expectCalled(want: bool, before: usize, what: []const u8) !void {
    if ((calls.load(.monotonic) > before) == want) return;
    check.say("{s}: Accelerate {s}called: ", .{ what, if (want) "not " else "" });
    return error.TestUnexpectedResult;
}

/// Counts calls into Accelerate routine by routine.
const Probe = struct {
    on: bool,
    last: usize,

    fn init(on_: bool) Probe {
        return .{ .on = on_, .last = calls.load(.monotonic) };
    }

    /// The routine just run called Accelerate if `routed` and `on`.
    fn saw(p: *Probe, routed: bool, what: []const u8) !void {
        try expectCalled(p.on and routed, p.last, what);
        p.last = calls.load(.monotonic);
    }
};

fn val(comptime T: type, re: Float(T), im: Float(T)) T {
    return if (comptime isComplex(T)) .init(re, im) else re;
}

/// Each public routine over `T` once, on small finite data: each one
/// Accelerate has calls it when `on`, and nothing else does.
fn everyRoutine(comptime T: type, on_: bool) !void {
    const al = zigstd.testAllocator;
    const io = check.io;
    const cx = comptime isComplex(T);
    const real = !cx;
    const al1 = val(T, 1.5, -0.5);
    const be1 = val(T, 0.5, 0.25);
    var a: [36]T = undefined;
    var bm: [36]T = undefined;
    var c: [36]T = undefined;
    var x: [6]T = undefined;
    var y: [6]T = undefined;
    check.fill(T, &a, 1);
    check.fill(T, &bm, 2);
    check.fill(T, &c, 3);
    check.fill(T, &x, 4);
    check.fill(T, &y, 5);
    // dominant diagonals, for the solves
    for (0..6) |i| a[i * 7] = val(T, 4, 0);
    for (0..6) |i| a[i] = val(T, 4, 0);
    const n: isize = 6;
    var p: Probe = .init(on_);

    try blas.copy(io, T, n, &x, 1, &y, 1);
    try p.saw(true, "copy");
    try blas.swap(io, T, n, &x, 1, &y, 1);
    try p.saw(true, "swap");
    // strided calls go by `sends`; see the test of it
    try blas.swap(io, T, n, &x, 1, &y, -1);
    try p.saw(x86_64, "swap, a negative stride");
    try blas.scal(io, T, n, al1, &x, 1);
    try p.saw(true, "scal");
    _ = try blas.dot(io, T, n, &x, 1, &y, 1);
    try p.saw(T != f32, "dot");
    if (cx) {
        _ = try blas.dotc(io, T, n, &x, 1, &y, 1);
        try p.saw(true, "dotc");
        try blas.axpyc(io, T, n, al1, &x, 1, &y, 1);
        try p.saw(false, "axpyc");
        // `alpha` has a zero imaginary part (`quirks`)
        try blas.scalReal(io, T, n, 0.5, &x, 1);
        try p.saw(false, "scalReal");
    } else {
        var param = [5]T{ -1, 0.5, 0.25, -0.25, 0.75 };
        try blas.rotm(io, T, n, &x, 1, &y, 1, &param);
        try p.saw(true, "rotm");
    }
    if (T == f32) {
        _ = try blas.dsdot(io, n, &x, 1, &y, 1);
        try p.saw(!x86_64, "dsdot");
        _ = try blas.sdsdot(io, n, 0.5, &x, 1, &y, 1);
        try p.saw(!x86_64, "sdsdot");
    }
    _ = try blas.nrm2(io, T, n, &x, 1);
    try p.saw(x86_64 and Float(T) == f64, "nrm2");
    _ = try blas.nrm2(io, T, n, &x, -1);
    try p.saw(false, "nrm2, a negative stride");
    _ = try blas.asum(io, T, n, &x, 1);
    try p.saw(true, "asum");
    _ = try blas.sum(io, T, n, &x, 1);
    try p.saw(false, "sum");
    _ = try blas.index(io, T, .amax, n, &x, 1);
    try p.saw(false, "index amax");
    _ = try blas.index(io, T, .amin, n, &x, 1);
    try p.saw(false, "index amin");
    _ = try blas.extreme(io, T, .amax, n, &x, 1);
    try p.saw(false, "extreme amax");
    try blas.rot(io, T, n, &x, 1, &y, 1, 0.6, 0.8);
    try p.saw(true, "rot");
    try blas.axpy(io, T, n, al1, &x, 1, &y, 1);
    try p.saw(true, "axpy");
    try blas.axpby(io, T, n, al1, &x, 1, be1, &y, 1);
    try p.saw(true, "axpby");

    try blas.gemv(al, io, T, .n, n, n, al1, &a, n, &x, 1, be1, &y, 1);
    try p.saw(true, "gemv");
    try blas.gemv(al, io, T, if (cx) .c else .t, n, n, al1, &a, n, &x, 1, be1, &y, 1);
    try p.saw(true, "gemv transposed");
    try blas.gemv(al, io, T, .r, n, n, al1, &a, n, &x, 1, be1, &y, 1);
    try p.saw(real, "gemv r");
    try blas.gbmv(al, io, T, .t, n, n, 1, 1, al1, &a, 3, &x, 1, be1, &y, 1);
    try p.saw(true, "gbmv");
    try blas.ger(al, io, T, n, n, al1, &x, 1, &y, 1, &c, n);
    try p.saw(true, "ger");
    if (cx) {
        try blas.gerc(al, io, T, n, n, al1, &x, 1, &y, 1, &c, n);
        try p.saw(true, "gerc");
        try blas.hemv(al, io, T, .l, n, al1, &a, n, &x, 1, be1, &y, 1);
        try p.saw(true, "hemv");
        try blas.hemv(al, io, T, .v, n, al1, &a, n, &x, 1, be1, &y, 1);
        try p.saw(false, "hemv v");
        try blas.hbmv(al, io, T, .u, n, 1, al1, &a, 2, &x, 1, be1, &y, 1);
        try p.saw(true, "hbmv");
        try blas.hpmv(al, io, T, .l, n, al1, &c, &x, 1, be1, &y, 1);
        try p.saw(true, "hpmv");
        try blas.her(al, io, T, .u, n, 0.5, &x, 1, &c, n);
        try p.saw(true, "her");
        try blas.hpr(al, io, T, .l, n, 0.5, &x, 1, &c);
        try p.saw(true, "hpr");
        try blas.her2(al, io, T, .u, n, al1, &x, 1, &y, 1, &c, n);
        try p.saw(true, "her2");
        try blas.hpr2(al, io, T, .l, n, al1, &x, 1, &y, 1, &c);
        try p.saw(true, "hpr2");
    }
    try blas.symv(al, io, T, .u, n, al1, &a, n, &x, 1, be1, &y, 1);
    try p.saw(real, "symv");
    try blas.sbmv(al, io, T, .l, n, 1, al1, &a, 2, &x, 1, be1, &y, 1);
    try p.saw(real, "sbmv");
    try blas.spmv(al, io, T, .u, n, al1, &c, &x, 1, be1, &y, 1);
    try p.saw(real, "spmv");
    try blas.syr(al, io, T, .u, n, al1, &x, 1, &c, n);
    try p.saw(real, "syr");
    try blas.spr(al, io, T, .l, n, al1, &x, 1, &c);
    try p.saw(real, "spr");
    try blas.syr2(al, io, T, .u, n, al1, &x, 1, &y, 1, &c, n);
    try p.saw(real, "syr2");
    try blas.spr2(al, io, T, .l, n, al1, &x, 1, &y, 1, &c);
    try p.saw(real, "spr2");
    try blas.trmv(al, io, T, .u, .t, .n, n, &a, n, &x, 1);
    try p.saw(true, "trmv");
    try blas.trmv(al, io, T, .u, .r, .n, n, &a, n, &x, 1);
    try p.saw(real, "trmv r");
    try blas.trsv(al, io, T, .l, .n, .u, n, &a, n, &x, -1);
    try p.saw(true, "trsv");
    try blas.tbmv(al, io, T, .u, .n, .n, n, 1, &a, 2, &x, 1);
    try p.saw(true, "tbmv");
    try blas.tbsv(al, io, T, .l, .c, .n, n, 1, &a, 2, &x, 1);
    try p.saw(true, "tbsv");
    try blas.tpmv(al, io, T, .u, .n, .u, n, &a, &x, 1);
    try p.saw(true, "tpmv");
    try blas.tpsv(al, io, T, .l, .n, .u, n, &a, &x, 1);
    try p.saw(true, "tpsv");

    try blas.gemm(al, io, T, .n, if (cx) .c else .t, n, n, n, al1, &a, n, &bm, n, be1, &c, n);
    try p.saw(true, "gemm");
    try blas.gemm(al, io, T, .r, .n, n, n, n, al1, &a, n, &bm, n, be1, &c, n);
    try p.saw(real, "gemm r");
    try blas.symm(al, io, T, .l, .u, n, n, al1, &a, n, &bm, n, be1, &c, n);
    try p.saw(true, "symm");
    try blas.syrk(al, io, T, .u, .n, n, n, al1, &a, n, be1, &c, n);
    try p.saw(true, "syrk");
    try blas.syr2k(al, io, T, .l, .t, n, n, al1, &a, n, &bm, n, be1, &c, n);
    try p.saw(true, "syr2k");
    if (cx) {
        try blas.hemm(al, io, T, .r, .l, n, n, al1, &a, n, &bm, n, be1, &c, n);
        try p.saw(true, "hemm");
        try blas.herk(al, io, T, .l, .c, n, n, 1.5, &a, n, 0.5, &c, n);
        try p.saw(true, "herk");
        try blas.her2k(al, io, T, .u, .n, n, n, al1, &a, n, &bm, n, 0.5, &c, n);
        try p.saw(true, "her2k");
    }
    try blas.trmm(al, io, T, .r, .u, .t, .n, n, n, al1, &a, n, &c, n);
    try p.saw(true, "trmm");
    try blas.trsm(al, io, T, .l, .l, .n, .u, n, n, al1, &a, n, &c, n);
    try p.saw(true, "trsm");
}

test "the public routines call Accelerate on macOS where it has them, for s, d, c and z only" {
    _ = blas.vtable();
    const was = routing();
    defer _ = route(was);
    if (!route(true)) return zigstd.testExpectTrue(!available or accel.loaded == null);
    inline for (.{ f32, f64, Complex32, Complex64 }) |T| try everyRoutine(T, true);
    // types it has nothing for
    var h: [4]f16 = @splat(1);
    const before = calls.load(.monotonic);
    _ = try blas.dot(check.io, f16, 4, &h, 1, &h, 1);
    inline for (check.x_real) |X| {
        var q: [4]X = @splat(1);
        _ = try blas.dot(check.io, X, 4, &q, 1, &q, 1);
    }
    try expectCalled(false, before, "f16 and XDouble dot");
}

test "the off switch runs rstd's kernels, and RSTD_BLAS_ACCELERATE=0 sets it at init" {
    _ = blas.vtable();
    const was = routing();
    defer _ = route(was);
    try zigstd.testExpectTrue(!blas.useAccelerate(false));
    inline for (.{ f32, f64, Complex32, Complex64 }) |T| try everyRoutine(T, false);
    if (!blas.useAccelerate(true)) return error.SkipZigTest;
    try everyRoutine(f64, true);

    blas.shutdown();
    defer {
        blas.shutdown();
        _ = blas.vtable();
    }
    var env: std.process.Environ.Map = .init(zigstd.testAllocator);
    defer env.deinit();
    try env.put("RSTD_BLAS_ACCELERATE", "0");
    try zigstd.testExpectTrue(blas.init(&env));
    try zigstd.testExpectTrue(!routing());
    try everyRoutine(Complex32, false);
    // loaded before, and never unloaded
    try zigstd.testExpectTrue(blas.useAccelerate(true));
    try everyRoutine(Complex32, true);
}

fn dotOnThread(out: *?Threading, before: usize, called: *bool) void {
    var x = [_]f64{ 1, 2, 3 };
    _ = blas.dot(check.io, f64, 3, &x, 1, &x, 1) catch {};
    called.* = calls.load(.monotonic) > before;
    out.* = accel.loaded.?.threading();
}

test "each thread that calls in runs Accelerate single-threaded" {
    if (sysinfo.single_threaded) return error.SkipZigTest;
    _ = blas.vtable();
    if (!routing()) return error.SkipZigTest;
    var model: ?Threading = .multi_threaded;
    var called = false;
    const th = try std.Thread.spawn(.{}, dotOnThread, .{ &model, calls.load(.monotonic), &called });
    th.join();
    try zigstd.testExpectTrue(called);
    // null before macOS 15, where only `VECLIB_MAXIMUM_THREADS` does it
    if (model) |m| try zigstd.testExpectEqual(Threading.single_threaded, m);
}

// ---- rstd's workers around Accelerate --------------------------------------------

/// A split test's sizes: under a `test_serial_max` of `split_unit`, several
/// parts each, every one with a tail.
const split_unit = vector_quantum;
const split_n = 5 * vector_quantum + 37;
const split_m = 3 * matrix_quantum + 11;
const split_k = 2 * matrix_quantum + 5;

/// A split test's outputs: their bits, and their values as `T`.
fn Outputs(comptime T: type) type {
    return struct {
        const Self = @This();
        bits: std.ArrayList(u8) = .empty,
        vals: std.ArrayList(T) = .empty,

        fn deinit(o: *Self, al: std.mem.Allocator) void {
            o.bits.deinit(al);
            o.vals.deinit(al);
        }

        fn add(o: *Self, al: std.mem.Allocator, v: []const T) !void {
            try o.bits.appendSlice(al, std.mem.sliceAsBytes(v));
            try o.vals.appendSlice(al, v);
        }

        /// A scalar result: `T`, or a real one of any float type.
        fn one(o: *Self, al: std.mem.Allocator, v: anytype) !void {
            try o.bits.appendSlice(al, std.mem.asBytes(&v));
            const V = @TypeOf(v);
            try o.vals.append(al, if (V == T) v else val(T, @floatCast(v), 0));
        }
    };
}

/// Every routine that splits, routed, over `T`, on `io`.
fn splitRun(comptime T: type, al: std.mem.Allocator, io: Io, o: *Outputs(T)) !void {
    const cx = comptime isComplex(T);
    const al1 = val(T, 1.25, -0.5);
    const be1 = val(T, 0.5, 0.25);
    const nl: isize = split_n;
    const x = try al.alloc(T, 2 * split_n);
    defer al.free(x);
    const y = try al.alloc(T, 3 * split_n);
    defer al.free(y);
    for ([_][2]isize{ .{ 1, 1 }, .{ -2, 3 } }) |inc| {
        check.fill(T, x, 1);
        check.fill(T, y, 2);
        const px: isize = @intCast(@abs(inc[0]));
        try blas.copy(io, T, nl, x.ptr, inc[0], y.ptr, inc[1]);
        try o.add(al, y);
        try blas.swap(io, T, nl, x.ptr, inc[0], y.ptr, inc[1]);
        try o.add(al, x);
        try blas.scal(io, T, nl, al1, x.ptr, px);
        try o.add(al, x);
        try blas.rot(io, T, nl, x.ptr, inc[0], y.ptr, inc[1], 0.6, -0.8);
        try o.add(al, x);
        try o.add(al, y);
        try blas.axpy(io, T, nl, al1, x.ptr, inc[0], y.ptr, inc[1]);
        try o.add(al, y);
        try o.one(al, try blas.dot(io, T, nl, x.ptr, inc[0], y.ptr, inc[1]));
        if (cx) try o.one(al, try blas.dotc(io, T, nl, x.ptr, inc[0], y.ptr, inc[1]));
        if (T == f32) try o.one(al, try blas.dsdot(io, nl, x.ptr, inc[0], y.ptr, inc[1]));
        try o.one(al, try blas.asum(io, T, nl, x.ptr, px));
    }

    const ld = split_m + 1;
    const len = ld * split_m;
    const a = try al.alloc(T, len);
    defer al.free(a);
    const b0 = try al.alloc(T, len);
    defer al.free(b0);
    const c = try al.alloc(T, len);
    defer al.free(c);
    check.fill(T, a, 3);
    check.fill(T, b0, 4);
    // a dominant diagonal, for the solves
    for (0..split_m) |i| a[i * (ld + 1)] = val(T, 2 * split_m, 0);
    const l: isize = ld;
    const big: isize = split_m;
    const small: isize = split_k;
    // the long side first, then the other way about: each side splits
    for ([_][2]isize{ .{ big, small }, .{ small, big } }) |mn| {
        const m = mn[0];
        const n = mn[1];
        for ([_]types.Trans{ .n, if (cx) .c else .t }) |t| {
            check.fill(T, y, 5);
            try blas.gemv(al, io, T, t, m, n, al1, a.ptr, l, x.ptr, -2, be1, y.ptr, 3);
            try o.add(al, y);
        }
        @memcpy(c, b0);
        try blas.ger(al, io, T, m, n, al1, x.ptr, 1, y.ptr, -3, c.ptr, l);
        try o.add(al, c);
        for ([_][2]types.Trans{ .{ .n, .n }, .{ .t, if (cx) .c else .t } }) |tt| {
            @memcpy(c, b0);
            try blas.gemm(al, io, T, tt[0], tt[1], m, n, small, al1, a.ptr, l, a.ptr + 1, l, be1, c.ptr, l);
            try o.add(al, c);
        }
        for ([_]types.Side{ .l, .r }) |sd| {
            @memcpy(c, b0);
            try blas.symm(al, io, T, sd, .u, m, n, al1, a.ptr, l, a.ptr + 1, l, be1, c.ptr, l);
            try o.add(al, c);
            if (cx) {
                @memcpy(c, b0);
                try blas.hemm(al, io, T, sd, .l, m, n, al1, a.ptr, l, a.ptr + 1, l, be1, c.ptr, l);
                try o.add(al, c);
            }
            @memcpy(c, b0);
            try blas.trmm(al, io, T, sd, .u, .t, .n, m, n, al1, a.ptr, l, c.ptr, l);
            try o.add(al, c);
            @memcpy(c, b0);
            try blas.trsm(al, io, T, sd, .l, .n, .n, m, n, al1, a.ptr, l, c.ptr, l);
            try o.add(al, c);
        }
    }
}

/// `splitRun`'s outputs at one worker against 3 and 8, on a threaded `Io`
/// and an inline one; then against unsplit calls, within rounding.
fn splitType(comptime T: type, al: std.mem.Allocator, tio: Io) !void {
    var want: Outputs(T) = .{};
    defer want.deinit(al);
    blas.config.cpu_count = 1;
    threading.test_serial_max = split_unit;
    const before = calls.load(.monotonic);
    try splitRun(T, al, Io.failing, &want);
    const split_calls = calls.load(.monotonic) - before;
    for ([_]u32{ 3, 8 }) |cpus| for ([_]Io{ tio, Io.failing }) |io| {
        blas.config.cpu_count = cpus;
        var got: Outputs(T) = .{};
        defer got.deinit(al);
        try splitRun(T, al, io, &got);
        try zigstd.testExpectEqualSlices(u8, want.bits.items, got.bits.items);
    };
    var whole_: Outputs(T) = .{};
    defer whole_.deinit(al);
    threading.test_serial_max = std.math.maxInt(u64);
    const before_whole = calls.load(.monotonic);
    try splitRun(T, al, Io.failing, &whole_);
    // the parts were Accelerate calls of their own
    try zigstd.testExpectTrue(split_calls > 2 * (calls.load(.monotonic) - before_whole));
    const tol = check.tolerance(Float(T), split_k);
    for (whole_.vals.items, want.vals.items, 0..) |w, g, i| {
        const scale = @max(1, @abs(if (comptime isComplex(T)) w.re else w) + @abs(if (comptime isComplex(T)) w.im else 0));
        check.expectClose(T, w, g, tol * scale) catch |e| {
            check.say("{s} output {d}: ", .{ @typeName(T), i });
            return e;
        };
    }
}

test "routed calls split by size on rstd's workers: one worker's bits at any worker count and Io" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    _ = blas.vtable();
    if (!routing()) return error.SkipZigTest;
    const cpus = blas.config.cpu_count;
    defer blas.config.cpu_count = cpus;
    defer threading.test_serial_max = null;
    const al = zigstd.testAllocator;
    var threaded: Io.Threaded = .init(al, .{});
    defer threaded.deinit();
    const before = unpinned.load(.monotonic);
    inline for (.{ f32, f64, Complex32, Complex64 }) |T| splitType(T, al, threaded.io()) catch |e| {
        check.say("{s}: ", .{@typeName(T)});
        return e;
    };
    // every worker ran its parts single-threaded (macOS 15+)
    try zigstd.testExpectEqual(before, unpinned.load(.monotonic));
}

test "on the 32-bit table, what does not fit c_int runs rstd's kernels" {
    // sizes past `c_int` need a 64-bit `isize`, as on every macOS
    if (comptime @bitSizeOf(isize) <= 32) return error.SkipZigTest;
    _ = blas.vtable();
    if (!routing()) return error.SkipZigTest;
    const lp = load(.{ .lp64 = true }) orelse return error.SkipZigTest;
    try zigstd.testExpectTrue(lp.blas == .lp64);
    const saved = accel.loaded;
    accel.loaded = lp;
    defer accel.loaded = saved;
    const io = check.io;
    const al = zigstd.testAllocator;
    const big: isize = std.math.maxInt(c_int) + 1;
    var x = [_]f64{ 1, -2, 3 };
    var y = [_]f64{ 0.5, 4, -1 };

    var before = calls.load(.monotonic);
    try zigstd.testExpectEqual(@as(f64, -10.5), try blas.dot(io, f64, 3, &x, 1, &y, 1));
    try expectCalled(true, before, "dot");
    try zigstd.testExpectEqual(@as(isize, 3), try blas.index(io, f64, .amax, 3, &x, 1));
    // one element at a stride past `c_int`
    before = calls.load(.monotonic);
    try zigstd.testExpectEqual(@as(f64, 0.5), try blas.dot(io, f64, 1, &x, big, &y, 1));
    try zigstd.testExpectEqual(@as(f64, 1), try blas.asum(io, f64, 1, &x, big));
    try zigstd.testExpectEqual(@as(isize, 1), try blas.index(io, f64, .amax, 1, &x, big));
    // one column at a leading dimension past it
    var c = [_]f64{ 1, 1, 1 };
    try blas.gemv(al, io, f64, .n, 3, 1, 2, &x, big, &y, 1, 1, &c, 1);
    try zigstd.testExpectEqualSlices(f64, &.{ 2, -1, 4 }, &c);
    try blas.gemm(al, io, f64, .n, .n, 3, 1, 1, 1, &x, 3, &y, big, 0, &c, big);
    try zigstd.testExpectEqualSlices(f64, &.{ 0.5, -1, 1.5 }, &c);
    try expectCalled(false, before, "past c_int");
    // an operand whose span is past `c_int`, though each argument fits
    const half: isize = std.math.maxInt(c_int) / 2 + 1;
    try zigstd.testExpectTrue(call("cblas_ddot", &accel.loaded.?.blas, &.{ span(3, half), 3 }, .{ @as(isize, 3), &x, half, &y, @as(isize, 1) }) == null);
    try zigstd.testExpectEqual(@as(f64, -10.5), call("cblas_ddot", &accel.loaded.?.blas, &.{ 3, 3 }, .{ @as(isize, 3), &x, @as(isize, 1), &y, @as(isize, 1) }).?);
}

// ---- Accelerate against rstd's kernels on special values ----------------------

/// The routed level-1 routines.
const R1 = enum { copy, swap, scal, dot, dotc, dsdot, nrm2, asum, iamax, rot, rotm, axpy, axpby };

/// `alpha`, `beta` (`scal`, `axpy`, `axpby`), `c`, `s` (`rot`) and `rotm`'s
/// flag, by variant.
const scalars1 = [_]struct { a: [2]f64, b: [2]f64, cs: [2]f64, flag: f64 }{
    .{ .a = .{ 1.5, -0.5 }, .b = .{ 0.5, 0.25 }, .cs = .{ 0.6, 0.8 }, .flag = -1 },
    .{ .a = .{ 0, 0 }, .b = .{ 0.5, 0 }, .cs = .{ 0, 1 }, .flag = 0 },
    .{ .a = .{ 1.5, 0 }, .b = .{ 0, 0 }, .cs = .{ 1, 0 }, .flag = 1 },
    .{ .a = .{ 0, 0 }, .b = .{ 0, 0 }, .cs = .{ -0.8, 0.6 }, .flag = -2 },
    .{ .a = .{ 1, 0 }, .b = .{ 1, 0 }, .cs = .{ 0.6, -0.8 }, .flag = -1 },
    .{ .a = .{ 0, 1 }, .b = .{ 0, -0.5 }, .cs = .{ 0.8, 0.6 }, .flag = 0 },
};

fn variants1(r: R1) usize {
    return switch (r) {
        .scal, .rot, .rotm, .axpy, .axpby => scalars1.len,
        else => 1,
    };
}

fn has1(comptime T: type, r: R1) bool {
    return switch (r) {
        .dotc => isComplex(T),
        .dsdot => T == f32,
        .rotm => !isComplex(T),
        else => true,
    };
}

/// A level-1 run's value or index.
const Out1 = struct { v: [2]f128 = .{ 0, 0 }, i: isize = 0 };

fn out1(v: anytype) Out1 {
    if (comptime isComplex(@TypeOf(v))) return .{ .v = .{ v.re, v.im } };
    return .{ .v = .{ v, 0 } };
}

fn run1(comptime T: type, comptime r: R1, k: usize, n: usize, x: *check.Vector(T), y: *check.Vector(T)) !Out1 {
    const io = check.io;
    const ni: isize = @intCast(n);
    const s = scalars1[k];
    const F = Float(T);
    const alpha = val(T, @floatCast(s.a[0]), @floatCast(s.a[1]));
    const beta = val(T, @floatCast(s.b[0]), @floatCast(s.b[1]));
    const xb = x.base();
    const yb = y.base();
    switch (r) {
        .copy => try blas.copy(io, T, ni, xb, x.inc, yb, y.inc),
        .swap => try blas.swap(io, T, ni, xb, x.inc, yb, y.inc),
        .scal => try blas.scal(io, T, ni, alpha, xb, x.inc),
        .dot => return out1(try blas.dot(io, T, ni, xb, x.inc, yb, y.inc)),
        .dotc => if (comptime isComplex(T)) return out1(try blas.dotc(io, T, ni, xb, x.inc, yb, y.inc)),
        .dsdot => if (comptime T == f32) return out1(try blas.dsdot(io, ni, xb, x.inc, yb, y.inc)),
        .nrm2 => return out1(try blas.nrm2(io, T, ni, xb, x.inc)),
        .asum => return out1(try blas.asum(io, T, ni, xb, x.inc)),
        .iamax => return .{ .i = try blas.index(io, T, .amax, ni, xb, x.inc) },
        .rot => try blas.rot(io, T, ni, xb, x.inc, yb, y.inc, @floatCast(s.cs[0]), @floatCast(s.cs[1])),
        .rotm => if (comptime !isComplex(T)) {
            const param = [5]F{ @floatCast(s.flag), 0.75, -0.5, 0.25, 1.25 };
            try blas.rotm(io, T, ni, xb, x.inc, yb, y.inc, &param);
        },
        .axpy => try blas.axpy(io, T, ni, alpha, xb, x.inc, yb, y.inc),
        .axpby => try blas.axpby(io, T, ni, alpha, xb, x.inc, beta, yb, y.inc),
    }
    return .{};
}

/// `a` and `b` alike: NaN for NaN, the same infinity, else within `ulps` of
/// `scale`; with `lax`, any non-finite for any (complex infinities meet in
/// orders of their own).
fn alike(comptime F: type, a: f128, b: f128, ulps: f128, scale: f128, lax: bool) bool {
    if (!std.math.isFinite(a) or !std.math.isFinite(b)) {
        if (lax and !std.math.isFinite(a) and !std.math.isFinite(b)) return true;
        if (std.math.isNan(a) or std.math.isNan(b)) return std.math.isNan(a) and std.math.isNan(b);
        return a == b;
    }
    return @abs(a - b) <= ulps * (std.math.floatEps(F) * scale + std.math.floatTrueMin(F));
}

/// Element `i` of `v` as the routine read it (a zero stride reads one).
fn input(comptime T: type, v: *const check.Vector(T), i: usize) T {
    return if (v.inc == 0) v.was[v.n - 1] else v.was[i];
}

/// The case's results against rstd's: values within their accuracy, an
/// index exactly, vectors element by element, gaps untouched.
fn compare1(comptime T: type, comptime r: R1, n: usize, want: Out1, got: Out1, x1: *const check.Vector(T), y1: *const check.Vector(T), x2: *const check.Vector(T), y2: *const check.Vector(T)) !void {
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
var mem1: [4 * (2 * 1025 * 3 + 64) * @sizeOf(Complex64)]u8 align(16) = undefined;

/// One case through rstd's kernels and through Accelerate, compared.
fn case1(comptime T: type, comptime r: R1, k: usize, n: usize, inc: [2]isize, p: ?pattern.Pattern, on_y: bool) !void {
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
const Seen = struct {
    keys: [256][96]u8 = undefined,
    lens: [256]usize = undefined,
    len: usize = 0,
    failed: bool = false,

    /// Records `key`; true the first time.
    fn first(s: *Seen, key: []const u8) bool {
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

fn patternName(p: ?pattern.Pattern) []const u8 {
    const q = p orelse return "plain";
    return switch (q) {
        .one => |o| @tagName(o.kind),
        else => @tagName(q),
    };
}

/// Every size, stride and pattern of `r` over `T`; a failure is reported once
/// per variant, pattern kind and stride sign.
fn sweep1(comptime T: type, comptime r: R1, seen: *Seen) void {
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
                    const key = std.fmt.bufPrint(&kb, "{s} {s} v{d} {s} in {s} inc {d},{d}", .{ @tagName(r), @typeName(T), k, patternName(p), if (on_y) "y" else "x", std.math.sign(inc[0]), std.math.sign(inc[1]) }) catch "?";
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
/// `sends` says: the `x86_64` `cblas_dsdot`, and the strides `nrm2`, `asum`
/// and `scal` (a quick return) leave out.
fn kept1(comptime T: type, r: R1, inc: isize) bool {
    return switch (r) {
        .dot, .dsdot => x86_64 and T == f32,
        .nrm2 => if (Float(T) == f32) x86_64 else inc <= 0,
        .asum, .scal => inc <= 0,
        .iamax => true,
        else => false,
    };
}

/// Runs `r` over `T` on 6 elements at stride `inc` in both vectors; fails
/// unless it called Accelerate exactly when `want`.
fn expect1(comptime T: type, comptime r: R1, inc: isize, want: bool) !void {
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
    try expect1(f32, .asum, 2, true);
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

    try blas.scal(io, Complex64, 4, .init(1.5, 0), &z, 1);
    try p.saw(false, "complex scal, alpha with a zero part");
    try zigstd.testExpectTrue(std.math.isNan(z[0].im));
    try blas.axpy(io, Complex64, 4, .init(0, 1), &z, 1, &w, 1);
    try p.saw(false, "complex axpy, alpha with a zero part");
    try blas.axpby(io, f64, 4, 0, &x, 1, 2, &y, 1);
    try p.saw(false, "axpby, zero alpha");
    try zigstd.testExpectEqualSlices(f64, &.{ 2, 0, 4, 1 }, &y);
    try blas.axpby(io, Complex64, 4, .init(1, 1), &z, 1, .init(2, 0), &w, 1);
    try p.saw(false, "complex axpby, beta with a zero part");
    try blas.trmv(al, io, f64, .u, .n, .n, 4, &a, 4, &x, 1);
    try p.saw(false, "trmv, x not finite");
    try blas.ger(al, io, f64, 4, 4, 1, &x, 1, &y, 1, &a, 4);
    try p.saw(false, "ger, x not finite");
    // a zero factor beside an infinity, through rstd; without one, Accelerate
    var sx = [_]f64{ 0, inf, -1, 2 };
    @memset(&a, 1);
    try blas.syr(al, io, f64, .l, 4, 1, &sx, 1, &a, 4);
    try p.saw(false, "syr, x(0) zero and x not finite");
    try zigstd.testExpectEqual(@as(f64, 1), a[1]);
    // `x(3)` and `y(3)` are zero, `x` at a negative stride
    var sy = [_]f64{ 1, 2, 3, -0.0 };
    try blas.syr2(al, io, f64, .u, 4, 1, &sx, -1, &sy, 1, &a, 4);
    try p.saw(false, "syr2, x(j) and y(j) zero and x not finite");
    sy[3] = 1;
    try blas.syr2(al, io, f64, .u, 4, 1, &sx, -1, &sy, 1, &a, 4);
    try p.saw(true, "syr2, no j with both zero");
    try blas.trmm(al, io, f64, .l, .l, .n, .n, 4, 4, 1, &a, 4, &b, 4);
    try p.saw(false, "trmm, B not finite");
    // the same cases, finite, call it
    x = .{ 1, 3, -1, 2 };
    @memset(&a, 1);
    try blas.trmv(al, io, f64, .u, .n, .n, 4, &a, 4, &x, 1);
    try p.saw(true, "trmv");
    @memset(&b, 0.5);
    try blas.trmm(al, io, f64, .l, .l, .n, .n, 4, 4, 1, &a, 4, &b, 4);
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
            try blas.her2k(al, io, Complex64, .u, .n, 3, 4, .init(0.5, 0.25), &h, 3, &h, 3, 0.5, &c, 3);
        } else {
            try blas.herk(al, io, Complex64, .l, .n, 3, 4, 1.5, &h, 3, 0.5, &c, 3);
        }
        try p.saw(true, "herk and her2k");
        for (0..3) |i| try zigstd.testExpectEqual(@as(f64, 0), c[i * 4].im);
    };
    // and `her2`, packed too
    var hx: [3]Complex64 = .{ .init(inf, 1), .init(1, 2), .init(-1, 0.5) };
    check.fill(Complex64, &c, 5);
    try blas.her2(al, io, Complex64, .l, 3, .init(0.5, 0.25), &hx, 1, &hx, 1, &c, 3);
    try p.saw(true, "her2");
    for (0..3) |i| try zigstd.testExpectEqual(@as(f64, 0), c[i * 4].im);
    check.fill(Complex64, &c, 5);
    try blas.hpr2(al, io, Complex64, .u, 3, .init(0.5, 0.25), &hx, 1, &hx, 1, &c);
    try p.saw(true, "hpr2");
    for ([_]usize{ 0, 2, 5 }) |i| try zigstd.testExpectEqual(@as(f64, 0), c[i].im);
    // `beta` 1 and nothing added leaves `C` as it was, through rstd
    check.fill(Complex64, &c, 4);
    for (0..3) |i| c[i * 4].im = 3;
    try blas.herk(al, io, Complex64, .l, .n, 3, 0, 1.5, &h, 3, 1, &c, 3);
    try p.saw(false, "herk, beta 1 and k 0");
    for (0..3) |i| try zigstd.testExpectEqual(@as(f64, 3), c[i * 4].im);
    // a non-finite `C` with a nonzero `beta`, through rstd; with a zero one, Accelerate
    check.fill(Complex64, &h, 3);
    check.fill(Complex64, &c, 4);
    c[1] = .init(1, inf);
    try blas.herk(al, io, Complex64, .l, .n, 3, 4, 1.5, &h, 3, 0.5, &c, 3);
    try p.saw(false, "herk, C not finite");
    try zigstd.testExpectTrue(!std.math.isNan(c[1].re));
    c[1] = .init(1, inf);
    try blas.herk(al, io, Complex64, .l, .n, 3, 4, 1.5, &h, 3, 0, &c, 3);
    try p.saw(true, "herk, beta 0");
}
