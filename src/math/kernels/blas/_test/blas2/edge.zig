//! Edge cases through every level-2 routine on every core this machine runs,
//! against an `f128` reference: `pattern.Pattern`'s special values one operand at
//! a time, NaN where a routine must not read; what may not propagate, `design.md`'s Verification.

const std = @import("std");
const blas = @import("../../__root__.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../check.zig");
const pattern = @import("../pattern.zig");
const portable = @import("../../portable.zig");
const accel = @import("../../apple_accel.zig");
const test_io = @import("../io.zig");
const cs = @import("edge/case.zig");
const rn = @import("edge/runs.zig");

const C = std.math.Complex;
const Float = check.FloatOf;
const isComplex = check.isComplex;
const Pattern = pattern.Pattern;
const Case = cs.Case;
const Routine = cs.Routine;
const Where = cs.Where;
const shapeOf = cs.shapeOf;
const trans_flags = cs.trans_flags;
const val = cs.val;
const mvRun = rn.mvRun;
const rankRun = rn.rankRun;
const trRun = rn.trRun;

/// A failure of the sweep, so one run names every failing routine; 0 for
/// none, else `@intFromError`. Atomic: the Accelerate sweep runs routines on
/// threads of their own.
pub var failed = std.atomic.Value(@Int(.unsigned, @bitSizeOf(anyerror))).init(0);

pub fn fail(e: anyerror) void {
    failed.store(@intFromError(e), .monotonic);
}

/// The failure since `failed` was cleared, if any.
pub fn failure() ?anyerror {
    const v = failed.load(.monotonic);
    return if (v == 0) null else @errorFromInt(v);
}

/// Runs `r` over its shapes, variants, operands and patterns, then its
/// special scalars, for `T`; names a failing case and goes on with the next
/// routine.
pub fn sweep(comptime T: type, comptime r: Routine) void {
    const complex = comptime isComplex(T);
    if (comptime r.complexOnly() and !complex) return;
    const gen = val(T, 1.25, -0.5);
    const beta = val(T, -0.75, 0.25);
    const nan = std.math.nan(Float(T));
    const zero = val(T, 0, 0);
    const wheres: []const Where = switch (comptime r.family()) {
        .mv => &.{ .a, .x, .y },
        .tr => &.{ .a, .x },
        .rank => if (r == .syr or r == .spr or r == .her or r == .hpr) &.{ .a, .x } else &.{ .a, .x, .y },
    };
    var turn: usize = 0;
    for (shapes()) |sh| for (0..r.variants(complex)) |v| {
        if (thinned(sh, v)) continue;
        const base = baseCase(r, sh);
        for (wheres) |where| {
            const s = shapeOf(r, v, base);
            const t = (r == .gemv or r == .gbmv) and trans_flags[v][0];
            const len = switch (where) {
                .a => s.entries(),
                .x => if (r.family() == .rank or t) s.m else s.n,
                .y => if (r.family() == .rank or t) s.n else s.m,
                .none => unreachable,
            };
            if (len == 0) continue;
            var buf: [13]Pattern = undefined;
            for (pattern.patternsFor(len, &buf)) |p| {
                // huge entries overflow a solve's partial sums, which the
                // reference's do not
                if (r.solve() and where == .a and p == .huge) continue;
                // nor do `x` times itself, where `x` holds huge values
                if ((r == .syr or r == .spr or r == .her or r == .hpr) and where == .x and (p == .huge or p == .mixed)) continue;
                var c = base;
                c.p = p;
                c.where = where;
                turn += 1;
                c.inc = if (turn % 2 == 0) .{ 1, 1 } else .{ -2, 3 };
                run(T, r, v, gen, beta, c) catch |e| {
                    report(e, T, r, v, c);
                    // Accelerate's run reports every failing case
                    if (!accel_run) return;
                };
            }
        }
        // the scalars that are not propagation, and NaN ones that are
        var c = base;
        c.inc = .{ -1, 2 };
        if (r.family() == .mv) for ([_][2]T{ .{ gen, zero }, .{ zero, beta }, .{ val(T, nan, 0), beta }, .{ gen, val(T, nan, 0) } }) |ab| {
            run(T, r, v, ab[0], ab[1], c) catch |e| {
                check.say("alpha {any}, beta {any}: ", .{ ab[0], ab[1] });
                return report(e, T, r, v, c);
            };
        };
        if (r.family() == .rank) for ([_]T{ zero, val(T, nan, 0) }) |a| {
            run(T, r, v, a, beta, c) catch |e| {
                check.say("alpha {any}: ", .{a});
                return report(e, T, r, v, c);
            };
        };
    };
    for (0..r.variants(complex)) |v| {
        var c = baseCase(r, large_shape);
        c.inc = large_incs[v % 2];
        c.pad = large_pad;
        run(T, r, v, gen, beta, c) catch |e| return report(e, T, r, v, c);
    }
}

/// m, n, kl, ku; the rigs take one, past the unrolls.
pub const all_shapes = [_][4]usize{ .{ 1, 1, 0, 0 }, .{ 7, 4, 1, 2 }, .{ 4, 17, 3, 0 }, .{ 33, 9, 2, 2 } };

/// The large strides' and leading dimension's shape (the rigs' one), strides
/// by variant, and padding: past a page of `f32`, and within the rigs' 32 MiB
/// heap at `Complex(f64)` (4.2 MB a vector).
pub const large_shape = all_shapes[2];
pub const large_incs = [_][2]isize{ .{ 4097, -4099 }, .{ -16411, 16417 } };
pub const large_pad = 4099;

/// Accelerate's: every shape, and one past its small-problem paths.
pub const accel_shapes = all_shapes ++ [_][4]usize{.{ 130, 67, 4, 6 }};
/// The largest order of any shape.
pub const max_order = 130;

/// Accelerate's largest shape takes every third variant: its `f128`
/// reference is slow.
pub fn thinned(sh: [4]usize, v: usize) bool {
    return sh[0] == max_order and v % 3 != 0;
}
pub var accel_run = false;

pub fn shapes() []const [4]usize {
    if (accel_run) return &accel_shapes;
    return if (check.rig()) all_shapes[2..3] else &all_shapes;
}

/// `r`'s case of shape `sh`, no pattern yet.
pub fn baseCase(comptime r: Routine, sh: [4]usize) Case {
    const general = r == .gemv or r == .gbmv or r == .ger or r == .gerc;
    // square routines take `n` from the shape's larger side
    const sq = @max(sh[0], sh[1]);
    return .{ .m = if (general) sh[0] else sq, .n = if (general) sh[1] else sq, .kl = sh[2], .ku = sh[3], .inc = .{ 1, 1 }, .p = .mixed, .where = .none };
}

/// Denormal results: denormal the operands that reach them (`x` tiny where
/// it meets itself, so its products are), the rest and the scalars near 1,
/// so a flushed input or result shows.
pub fn denormals(comptime T: type, comptime r: Routine) void {
    const complex = comptime isComplex(T);
    if (comptime r.complexOnly() and !complex) return;
    const self = r == .syr or r == .spr or r == .her or r == .hpr;
    for (shapes()) |sh| for (0..r.variants(complex)) |v| for ([_][2]isize{ .{ 1, 1 }, .{ -2, 3 } }) |inc| {
        if (thinned(sh, v)) continue;
        var c = baseCase(r, sh);
        c.inc = inc;
        c.p = .denormal;
        switch (comptime r.family()) {
            .mv => {
                c.where = .a;
                c.also = .{ .where = .y, .p = .denormal };
            },
            .tr => c.where = .x,
            .rank => {
                c.where = .a;
                c.also = .{ .where = .x, .p = if (self) .tiny else .denormal };
            },
        }
        run(T, r, v, val(T, 1.25, -0.5), val(T, -0.75, 0.25), c) catch |e| return report(e, T, r, v, c);
    };
}

pub fn run(comptime T: type, comptime r: Routine, v: usize, alpha: T, beta: T, c: Case) !void {
    return switch (comptime r.family()) {
        .mv => mvRun(T, r, v, alpha, beta, c),
        .tr => trRun(T, r, v, c),
        .rank => rankRun(T, r, v, alpha, c),
    };
}

pub fn report(e: anyerror, comptime T: type, comptime r: Routine, v: usize, c: Case) void {
    if (c.where != .none) c.p.describe();
    check.say("in {s}, {s} {s} variant {d} m={d} n={d} k={d},{d} inc={d},{d} pad={d}\n", .{ @tagName(c.where), @tagName(r), @typeName(T), v, c.m, c.n, c.kl, c.ku, c.inc[0], c.inc[1], c.pad });
    fail(e);
}

/// Runs `jobs` on up to a thread per CPU, each taking the next until none
/// remain; in turn on this thread where there are no threads.
pub fn onThreads(comptime jobs: []const *const fn () void) void {
    const Queue = struct {
        next: std.atomic.Value(usize) = .init(0),

        pub fn work(q: *@This()) void {
            while (true) {
                const i = q.next.fetchAdd(1, .monotonic);
                if (i >= jobs.len) return;
                jobs[i]();
            }
        }
    };
    var q: Queue = .{};
    if (comptime !test_io.hosted) {
        q.work();
    } else {
        var threads: [jobs.len]std.Thread = undefined;
        const cpus = std.Thread.getCpuCount() catch 1;
        var started: usize = 0;
        defer for (threads[0..started]) |t| t.join();
        while (started + 1 < @min(cpus, jobs.len)) : (started += 1) {
            threads[started] = std.Thread.spawn(.{}, Queue.work, .{&q}) catch break;
        }
        q.work();
    }
}

/// Real and complex apart, and the complex ones on every core in two (the
/// products, then the rest): one test of more is too big an image for
/// sparc64's debug rig.
pub const real_types = .{ f32, f64 };
pub const complex_types = .{ C(f32), C(f64) };
pub const all_routines = std.enums.values(Routine);
pub const products = routinesOf(true);
pub const solves_and_updates = routinesOf(false);

/// The routines that are products (`mv`), or the rest.
pub fn routinesOf(comptime mv: bool) []const Routine {
    comptime {
        var n = 0;
        for (all_routines) |r| n += @intFromBool((r.family() == .mv) == mv);
        var out: [n]Routine = undefined;
        var i = 0;
        for (all_routines) |r| if ((r.family() == .mv) == mv) {
            out[i] = r;
            i += 1;
        };
        const final = out;
        return &final;
    }
}

/// `vt`'s kernels, with its own `dtb_entries` and with 3, which takes the
/// triangular drivers through `gemv`; static, as a vtable is large for a
/// rig's stack.
pub var small_dtb: Vtable = undefined;

pub fn CheckAll(comptime types: anytype, comptime routines: []const Routine) type {
    return struct {
        pub fn sweepAll(vt: *const Vtable) !void {
            failed.store(0, .monotonic);
            {
                const saved = blas.testUse(vt);
                defer blas.testRestore(saved);
                inline for (types) |T| inline for (routines) |r| sweep(T, r);
            }
            small_dtb = vt.*;
            small_dtb.dtb_entries = 3;
            const saved = blas.testUse(&small_dtb);
            defer blas.testRestore(saved);
            inline for (types) |T| inline for (routines) |r| if (comptime r == .trmv or r == .trsv) sweep(T, r);
            if (failure()) |e| return e;
        }

        pub fn denormalsAll(vt: *const Vtable) !void {
            failed.store(0, .monotonic);
            const saved = blas.testUse(vt);
            defer blas.testRestore(saved);
            inline for (types) |T| inline for (routines) |r| denormals(T, r);
            if (failure()) |e| return e;
        }

        /// Each type's routine, swept and with denormals.
        const jobs = blk: {
            var j: []const *const fn () void = &.{};
            for (types) |T| for (routines) |r| {
                if (r.complexOnly() and !isComplex(T)) continue;
                j = j ++ .{&struct {
                    fn run() void {
                        sweep(T, r);
                        denormals(T, r);
                    }
                }.run};
            };
            break :blk j;
        };

        /// Through Apple Accelerate, where the routines go to it.
        pub fn accelerate() !void {
            if (comptime !accel.covers(f64)) return error.SkipZigTest;
            const saved = blas.testUse(blas.vtable());
            defer blas.testRestore(saved);
            if (!accel.route(true)) return error.SkipZigTest;
            return accelSweep();
        }

        /// `accelerate`'s sweep, on whatever the table routes to.
        pub fn accelSweep() !void {
            failed.store(0, .monotonic);
            accel_run = true;
            defer accel_run = false;
            onThreads(jobs);
            if (failure()) |e| return e;
        }

        pub fn portableKernels() !void {
            var vt: Vtable = .{};
            portable.install(&vt);
            try sweepAll(&vt);
            try denormalsAll(&vt);
        }
    };
}

test "portable level-2 routines, real: special values, and NaN wherever they must not read or write" {
    try CheckAll(real_types, all_routines).portableKernels();
}

test "portable level-2 routines, complex: special values, and NaN wherever they must not read or write" {
    try CheckAll(complex_types, all_routines).portableKernels();
}

test "every linked core's level-2 routines, real: special values, and NaN wherever they must not read or write" {
    try check.eachCore(CheckAll(real_types, all_routines).sweepAll);
}

test "every linked core's level-2 routines, complex products: special values, and NaN wherever they must not read or write" {
    try check.eachCore(CheckAll(complex_types, products).sweepAll);
}

test "every linked core's level-2 routines, complex solves and updates: special values, and NaN wherever they must not read or write" {
    try check.eachCore(CheckAll(complex_types, solves_and_updates).sweepAll);
}

test "the level-2 routines through Accelerate, real: special values, denormals, and NaN wherever they must not read or write" {
    try CheckAll(real_types, all_routines).accelerate();
}

test "the level-2 routines through Accelerate, complex: special values, denormals, and NaN wherever they must not read or write" {
    try CheckAll(complex_types, all_routines).accelerate();
}

test "every linked core's level-2 routines, real: denormal results" {
    try check.eachCore(CheckAll(real_types, all_routines).denormalsAll);
}

test "every linked core's level-2 routines, complex products: denormal results" {
    try check.eachCore(CheckAll(complex_types, products).denormalsAll);
}

test "every linked core's level-2 routines, complex solves and updates: denormal results" {
    try check.eachCore(CheckAll(complex_types, solves_and_updates).denormalsAll);
}
