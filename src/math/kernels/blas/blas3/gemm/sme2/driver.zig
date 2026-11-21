//! The SME2 direct real `gemm` of Apple M4 and later (`vortexm4`): there
//! upstream's `interface/gemm.c` (clang builds) sends real `gemm` to
//! `SME_SGEMM_KERNEL` and `SME_DGEMM_KERNEL`, whose `s2_gemm` and `s2_run`
//! (OpenBLAS `kernel/arm64/sme2_gemm_impl.h`, `v0.3.34-794-gdbf2c24c9`,
//! BSD-3-Clause, the OpenBLAS Project) are ported here: plain loops for small
//! jobs, the transposed core problem, its parts, their blocking and the
//! workspace, from the caller's allocator where upstream skips the multiply
//! when its buffers fail (#593). The `vortexm4` object's slot (`kernel.zig`)
//! runs one part. Each element of `C` is one fused chain over the depth from
//! `beta C` (0 for a zero `beta`), `alpha` folded into the packed `op(B)`,
//! whatever the blocking and parts: any worker count gives the same bits.

const std = @import("std");
const sysinfo = @import("builtin");
const Vtable = @import("../../../Vtable.zig");
const threading = @import("../../../threading.zig");
const osinfo = @import("../../../cpu/osinfo.zig");
const job = @import("job.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;

pub const Error = Allocator.Error || Io.Cancelable;

/// Plain loops below this many multiply-adds, where streaming mode costs more
/// (`S2_SMALL_WORK`, Apple M4, both precisions).
const small_work = 3000;
/// Parts from this many multiply-adds (`1 << 22`).
const parts_work = 1 << 22;
/// The most parts (`S2_MAX_PARTS`).
const parts_max = 4;
/// The kernels keep sizes under this (`S2_FITS`).
const fits = 1 << 30;

/// `C := alpha*op(A)*op(B) + beta*C` for `f32` and `f64` where `vt` has the
/// SME2 kernel and this CPU runs it; false where not, for the level-3 driver.
/// `k` and `alpha` are nonzero (`driver.gemm`); up to `want` workers.
pub fn gemm(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, ta: bool, tb: bool, m: usize, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize, want: u32) Error!bool {
    if (comptime T != f32 and T != f64) return false;
    const kernel = (if (T == f32) vt.sme_sgemm_kernel else vt.sme_dgemm_kernel) orelse return false;
    if (m >= fits or n >= fits or k >= fits or !usable()) return false;
    const work = @as(f64, @floatFromInt(m)) * @as(f64, @floatFromInt(n)) * @as(f64, @floatFromInt(k));
    if (work < small_work) {
        small(T, ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        return true;
    }
    // column-major `C` as the row-major `C^T = op(B)^T op(A)^T`
    const whole: job.Job(T) = .{
        .m = n,
        .n = m,
        .k = k,
        .alpha = alpha,
        .beta = beta,
        .a = b,
        .lda = ldb,
        .a_cols = tb,
        .b = a,
        .ldb = lda,
        .b_cols = ta,
        .c = c,
        .ldc = ldc,
        .blk = undefined,
        .ac = undefined,
        .bc = undefined,
    };
    var jobs: [parts_max]job.Job(T) = undefined;
    const used = split(T, whole, if (work >= parts_work) partCount(T, whole, want) else 1, &jobs);
    var at: [parts_max + 1]usize = undefined;
    at[0] = 0;
    for (jobs[0..used], 0..) |*p, i| {
        p.blk = job.blocking(T, p.m, p.n, p.k);
        at[i + 1] = at[i] + job.Work.of(T, p.blk).len;
    }
    const buf = try allocator.alignedAlloc(T, .fromByteUnits(128), at[used]);
    defer allocator.free(buf);
    for (jobs[0..used], 0..) |*p, i| {
        p.ac = buf.ptr + at[i];
        p.bc = p.ac + job.Work.of(T, p.blk).a_len;
    }
    if (used == 1) {
        kernel(&jobs[0]);
        return true;
    }
    var plan: threading.Plan = .{ .chunks = .{ .count = used, .bounds = undefined }, .workers = @intCast(used) };
    for (0..used + 1) |i| plan.chunks.bounds[i] = i;
    const Parts = struct {
        kernel: @TypeOf(kernel),
        jobs: []const job.Job(T),

        fn one(s: @This(), ch: threading.Chunk) void {
            for (s.jobs[ch.start..ch.end]) |*p| s.kernel(p);
        }
    };
    try threading.forEach(io, &plan, Parts{ .kernel = kernel, .jobs = jobs[0..used] }, Parts.one);
    return true;
}

/// One part per SME unit, as many as `want` workers, each at least four
/// panels wide along the split (`s2_gemm`).
fn partCount(comptime T: type, jb: job.Job(T), want: u32) usize {
    const big = @max(jb.m, jb.n);
    return std.math.clamp(@min(units(), want, big / (4 * job.Shape(T).vl), parts_max), 1, parts_max);
}

/// `jb` cut into `parts` along its longer side, on whole groups of four
/// panels (`s2_gemm`); returns how many.
fn split(comptime T: type, jb: job.Job(T), parts: usize, out: *[parts_max]job.Job(T)) usize {
    const split_m = jb.m >= jb.n;
    const total = if (split_m) jb.m else jb.n;
    const quantum = 4 * job.Shape(T).vl;
    const step = (((total + parts - 1) / parts) + quantum - 1) / quantum * quantum;
    var used: usize = 0;
    var off: usize = 0;
    while (used < parts and off < total) : (off += step) {
        const len = @min(total - off, step);
        var p = jb;
        if (split_m) {
            p.m = len;
            p.a = if (jb.a_cols) jb.a + off else jb.a + off * jb.lda;
            p.c = jb.c + off * jb.ldc;
        } else {
            p.n = len;
            p.b = if (jb.b_cols) jb.b + off * jb.ldb else jb.b + off;
            p.c = jb.c + off;
        }
        out[used] = p;
        used += 1;
    }
    return used;
}

/// Column-major `gemm` in plain loops (`s2_small`): `beta` first, then each
/// column's sums, fused as clang contracts upstream's.
fn small(comptime T: type, ta: bool, tb: bool, m: usize, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize) void {
    // `B`'s stride along the depth and along the columns
    const sbl = if (tb) ldb else 1;
    const sbj = if (tb) 1 else ldb;
    for (0..n) |j| {
        const cj = c + j * ldc;
        if (beta == 0) @memset(cj[0..m], 0) else if (beta != 1) for (cj[0..m]) |*x| {
            x.* *= beta;
        };
    }
    if (ta) {
        for (0..n) |j| for (0..m) |i| {
            const ai = a + i * lda;
            const bj = b + j * sbj;
            var s: T = 0;
            for (0..k) |l| s = @mulAdd(T, ai[l], bj[l * sbl], s);
            c[i + j * ldc] = @mulAdd(T, alpha, s, c[i + j * ldc]);
        };
        return;
    }
    var j: usize = 0;
    while (j + 4 <= n) : (j += 4) {
        const c0 = c + j * ldc;
        const c1 = c0 + ldc;
        const c2 = c1 + ldc;
        const c3 = c2 + ldc;
        const bj = b + j * sbj;
        for (0..k) |l| {
            const al = a + l * lda;
            const bl = bj + l * sbl;
            const t0 = alpha * bl[0];
            const t1 = alpha * bl[sbj];
            const t2 = alpha * bl[2 * sbj];
            const t3 = alpha * bl[3 * sbj];
            for (0..m) |i| {
                const x = al[i];
                c0[i] = @mulAdd(T, t0, x, c0[i]);
                c1[i] = @mulAdd(T, t1, x, c1[i]);
                c2[i] = @mulAdd(T, t2, x, c2[i]);
                c3[i] = @mulAdd(T, t3, x, c3[i]);
            }
        }
    }
    while (j < n) : (j += 1) {
        const cj = c + j * ldc;
        for (0..k) |l| {
            const t = alpha * b[j * sbj + l * sbl];
            const al = a + l * lda;
            for (0..m) |i| cj[i] = @mulAdd(T, t, al[i], cj[i]);
        }
    }
}

/// 0 unknown, 1 no, 2 yes.
var usable_state: std.atomic.Value(u8) = .init(0);

/// Whether this CPU runs the kernels: FEAT_SME2 and a 512-bit streaming
/// vector (`s2_usable`), asked of the OS once.
pub fn usable() bool {
    return switch (usable_state.load(.monotonic)) {
        1 => false,
        2 => true,
        else => blk: {
            const ok = detect();
            usable_state.store(if (ok) 2 else 1, .monotonic);
            break :blk ok;
        },
    };
}

fn detect() bool {
    if (comptime sysinfo.target.cpu.arch != .aarch64) return false;
    if (comptime sysinfo.target.os.tag.isDarwin()) {
        return (osinfo.sysctlInt("hw.optional.arm.FEAT_SME2") orelse 0) != 0 and osinfo.sysctlInt("hw.optional.arm.sme_max_svl_b") == job.svl_bytes;
    }
    if (comptime sysinfo.target.os.tag == .linux) {
        const hwcap2_sme2: usize = 1 << 37;
        if (osinfo.hwcap2() & hwcap2_sme2 == 0) return false;
        // `PR_SME_GET_VL`: the streaming vector length in its low 16 bits
        const r = std.os.linux.prctl(64, 0, 0, 0, 0);
        return r < 1 << 32 and r & 0xffff == job.svl_bytes;
    }
    return false;
}

/// Tests set this to stand in for `units`.
pub var test_units: ?usize = null;

/// SME units threads can use (`s2_units`): Apple M4 and M5 have one per
/// cluster of the fastest cores; with a second tier named "Performance"
/// (M5 Pro, Max, the M6 here), two per such cluster, upstream's reading of
/// a third-party measurement. One elsewhere.
fn units() usize {
    if (sysinfo.is_test) if (test_units) |u| return u;
    var u: usize = 1;
    if (comptime sysinfo.target.os.tag.isDarwin() and sysinfo.link_libc) {
        const cpus = osinfo.sysctlInt("hw.perflevel0.physicalcpu") orelse 0;
        const per = osinfo.sysctlInt("hw.perflevel0.cpusperl2") orelse 0;
        if (per > 0 and cpus / per > 1) u = @intCast(cpus / per);
        var name: [32]u8 = @splat(0);
        var len: usize = name.len - 1;
        if (std.c.sysctlbyname("hw.perflevel1.name", &name, &len, null, 0) == 0 and std.mem.eql(u8, std.mem.sliceTo(&name, 0), "Performance")) u *= 2;
    }
    return u;
}
