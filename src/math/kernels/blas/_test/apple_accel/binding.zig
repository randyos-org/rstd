//! Binding: the tables each macOS binds, agreement with rstd, and reproducibility on one thread.

const std = @import("std");
const sysinfo = @import("builtin");
const dll = @import("../../../../../dll.zig");
const types = @import("../../types.zig");
const blas = @import("../../__root__.zig");
const l3 = @import("../../blas3/core.zig");
const l1 = @import("../../blas1/interface.zig");
const accel = @import("../../apple_accel.zig");
const zigstd = @import("../../../../../zigstd.zig");
const check = @import("../check.zig");
const util = @import("util.zig");

const Complex64 = types.Complex64;
const AppleAccelerate = accel.AppleAccelerate;
const Blas = accel.Blas;
const Table = accel.Table;
const ThreadApi = accel.ThreadApi;
const Threading = accel.Threading;
const available = accel.available;
const bindBlas = accel.bindBlas;
const fortran = accel.fortran;
const ilp64Names = accel.ilp64Names;
const load = accel.load;
const lp64Names = accel.lp64Names;
const route = accel.route;
const routing = accel.routing;
const fill = util.fill;

/// The integer type of `tag`'s table.
pub fn Int(comptime tag: std.meta.Tag(Blas)) type {
    return if (tag == .ilp64) i64 else c_int;
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
/// LAPACK one (`?potrf_`, `?getrf_`, `?getrs_`) also without its underscore;
/// 13.3 adds each with `$NEWLAPACK` and
/// `$NEWLAPACK$ILP64`, LAPACK's without the underscore. None has
/// `BLASSetThreading` (macOS 15+). Re-derive from a newer SDK's `.tbd`s
/// rather than editing by hand.
pub const RecordedAccelerate = struct {
    /// macOS 13.3 or later
    new: bool,

    pub fn has(r: RecordedAccelerate, name: []const u8) bool {
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
        // 166 fields, less the four `i?amin`
        try zigstd.testExpectEqual(@as(usize, 162), b.lp64.count());
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
    try zigstd.testExpectEqual(@as(usize, 162), ilp.ilp64.count());
    try zigstd.testExpectEqualStrings("cblas_dgemm$NEWLAPACK$ILP64", ilp.ilp64.name(ilp64Names, .cblas_dgemm).?);
    try zigstd.testExpectEqualStrings("catlas_daxpby$NEWLAPACK$ILP64", ilp.ilp64.name(ilp64Names, .catlas_daxpby).?);
    try zigstd.testExpectEqualStrings("dpotrf$NEWLAPACK$ILP64", ilp.ilp64.name(ilp64Names, .dpotrf_).?);
    const lp = bindBlas(&new, true) orelse return error.TestUnexpectedResult;
    try zigstd.testExpectTrue(lp == .lp64);
    try zigstd.testExpectEqual(@as(usize, 162), lp.lp64.count());
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
            const want = try l1.dot(check.io, f64, 300, &x, 2, &y, -1);
            try zigstd.testExpectApproxEqRel(want, got, 1e-12);

            // zdotu_sub and zdotc_sub
            var zx: [257]Complex64 = undefined;
            var zy: [257]Complex64 = undefined;
            fill(f64, rng, @as(*[514]f64, @ptrCast(&zx)));
            fill(f64, rng, @as(*[514]f64, @ptrCast(&zy)));
            var zd: Complex64 = undefined;
            t.cblas_zdotu_sub.?(zx.len, &zx, 1, &zy, 1, &zd);
            const zw = try l1.dot(check.io, Complex64, zx.len, &zx, 1, &zy, 1);
            try zigstd.testExpectApproxEqAbs(zw.re, zd.re, 1e-12);
            try zigstd.testExpectApproxEqAbs(zw.im, zd.im, 1e-12);
            t.cblas_zdotc_sub.?(zx.len, &zx, 1, &zy, 1, &zd);
            const zc = try l1.dotc(check.io, Complex64, zx.len, &zx, 1, &zy, 1);
            try zigstd.testExpectApproxEqAbs(zc.re, zd.re, 1e-12);
            try zigstd.testExpectApproxEqAbs(zc.im, zd.im, 1e-12);

            // isamax: 0-based where rstd's is 1-based
            var sx: [1000]f32 = undefined;
            fill(f32, rng, &sx);
            sx[613] = -3;
            const ia: I = t.cblas_isamax.?(sx.len, &sx, 1);
            try zigstd.testExpectEqual(@as(I, 613), ia);
            try zigstd.testExpectEqual(try l1.index(check.io, f32, .largest_magnitude, sx.len, &sx, 1), ia + 1);

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
            try l3.gemm(gpa, check.io, f64, .t, .n, m, n, k, 1.5, a.ptr, k, bm.ptr, k, -0.5, c1.ptr, m);
            for (c0, c1) |g, w| try zigstd.testExpectApproxEqAbs(w, g, 1e-12 * k);
        },
    }
}

/// Each routine's result bits.
pub const Digest = struct { ddot: u64, dnrm2: u64, dasum: u64, dgemv_n: u64, dgemv_t: u64, dgemm: u64 };

pub const ReproData = struct {
    x: []f64,
    y: []f64,
    a: []f64,
    b: []f64,
    c: []f64,
    v: []f64,
    const n1 = 1_000_000;
    const n2 = 2000;
};

pub fn hash(s: []const f64) u64 {
    return std.hash.Wyhash.hash(0, std.mem.sliceAsBytes(s));
}

pub fn digest(comptime I: type, t: *const Table(I), d: *const ReproData) Digest {
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

pub fn digestOnThread(comptime I: type, acc: *const AppleAccelerate, t: *const Table(I), d: *const ReproData, out: *Digest) void {
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
