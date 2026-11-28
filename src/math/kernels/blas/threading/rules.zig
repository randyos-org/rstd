//! The thresholds (`threading.zig`): each routine family's rule, OpenBLAS's.

const std = @import("std");
const types = @import("../types.zig");
const coretypes = @import("../cpu/coretypes.zig");

const chunks_max = @import("chunks.zig").chunks_max;

/// OpenBLAS's `GEMM_MULTITHREAD_THRESHOLD` default.
const mt = 4;

/// A routine family, as OpenBLAS thresholds it, with the `work` its rule reads.
/// OpenBLAS runs a zero-stride level-1 call serially; so must the routine.
pub const Family = enum {
    /// `axpy`, `axpby`: `n`.
    axpy,
    /// `scal`: `n`.
    scal,
    /// `swap`: `n`.
    swap,
    /// `copy`: `n`.
    copy,
    /// `rot`, `rotm`: `n`.
    rot,
    /// `dot`, `dotc` and the mixed dots: `n`.
    dot,
    /// `asum`, `sum`: `n`.
    asum,
    /// `index`, `extreme`: `n`.
    index,
    /// `gemv`: `m * n`.
    gemv,
    /// `gbmv`: `m * n`.
    gbmv,
    /// `ger`: `m * n`.
    ger,
    /// `symv`, `spmv`, `hpmv`: `n`.
    symv,
    /// `hemv`: `n`.
    hemv,
    /// `trmv`: `n`.
    trmv,
    /// `tbmv`, `tpmv`: `n`.
    tpmv,
    /// `syr`, `spr`, `syr2`, `spr2`, real and complex symmetric: `n`.
    syr,
    /// `her`, `hpr`, `her2`, `hpr2`: `n`.
    her,
    /// `gemm`: `m * n * k`.
    gemm,
    /// `syrk`, `herk`: `(n + 1) * n * k`.
    syrk,
    /// `syr2k`, `her2k`: `n * k`.
    syr2k,
    /// `symm`: `2 * m * m * n`.
    symm,
    /// `trsm`: `m * n`.
    trsm,
    /// `omatcopy`, real: `rows * cols`.
    omatcopy,
};

/// Work below `below` takes at most `max` workers.
pub const Tier = struct { below: u64, max: u32 };

/// Work up to `x` takes at most `max` workers.
fn upTo(x: u64, max: u32) Tier {
    return .{ .below = x + 1, .max = max };
}

/// A family's threshold for one element type and core.
pub const Rule = struct {
    /// The most work run serially, as one chunk.
    serial_max: u64,
    /// Past the tiers, at most `work / serial_max` workers rather than any.
    per_unit: bool = false,
    /// Caps above `serial_max`, by ascending `below`.
    tiers: []const Tier = &.{},
    /// Units a chunk past `serial_max`; 0 for `serial_max`.
    chunk: u64 = 0,

    /// The most workers `work` justifies.
    pub fn sizeMax(r: Rule, work: u64) u32 {
        if (work <= r.serial_max) return 1;
        for (r.tiers) |t| if (work < t.below) return t.max;
        if (!r.per_unit) return std.math.maxInt(u32);
        return std.math.cast(u32, work / @max(r.serial_max, 1)) orelse std.math.maxInt(u32);
    }

    /// One chunk up to `serial_max`; above, chunks of about `chunk` (else
    /// `serial_max`), at most `chunks_max`.
    pub fn chunkCount(r: Rule, work: u64) usize {
        if (work <= r.serial_max) return 1;
        const per = @max(if (r.chunk != 0) r.chunk else r.serial_max, 1);
        return @intCast(@min(work / per + @intFromBool(work % per != 0), chunks_max));
    }
};

/// rstd's own level-2 thresholds on `haswell` (Coffee Lake, the Intel Mac), a
/// worker per call from each order up: where a second one paid in the bench's
/// `l2thread` runs (issue #845), against OpenBLAS's 50 to 100, below which
/// rstd's pool costs more than the call. The complex symmetric `syr` keeps
/// OpenBLAS's, untimed.
fn haswellLevel2(comptime family: Family, comptime T: type, core: anytype) ?Rule {
    if (!isX86(core, .haswell)) return null;
    const complex = comptime types.isComplex(T);
    const F = types.Float(T);
    return switch (family) {
        .syr => if (complex) null else .{ .serial_max = 511 },
        .tpmv => .{ .serial_max = 511 },
        .her => .{ .serial_max = 255 },
        .trmv => .{ .serial_max = if (!complex) 511 else if (F == f32) 383 else 255 },
        else => null,
    };
}

/// `family`'s rule for `T` on `core` (`blas.coreType()`; any architecture's
/// core enum, or null for none): OpenBLAS's
/// interface or kernel thresholds, with its per-core overrides as a
/// `DYNAMIC_ARCH` build applies them.
pub fn rule(comptime family: Family, comptime T: type, core: anytype) Rule {
    if (haswellLevel2(family, T, core)) |r| return r;
    const complex = comptime types.isComplex(T);
    const F = types.Float(T);
    const s: u64 = @sizeOf(F);
    return switch (family) {
        .axpy => .{ .serial_max = 10000 },
        .scal => .{ .serial_max = 1048576 },
        .swap => .{ .serial_max = @as(u64, if (complex) 1048576 else 2097152) * mt / s - 1 },
        // `kernel/arm64/copy_thunderx2t99.c`: `n <= 10000`
        .copy => .{ .serial_max = 10000 },
        .rot => .{ .serial_max = 100000 },
        .dot => (if (complex) null else dotCore(F, core)) orelse .{ .serial_max = 10000 },
        .asum => asumCore(core) orelse .{ .serial_max = if (complex) 10000 else 100000, .per_unit = true },
        // `kernel/arm64/iamax_thunderx2t99.c`, `izamax_thunderx2t99.c`: `n <= 10000`
        .index => .{ .serial_max = 10000 },
        .gemv => if (complex) .{ .serial_max = 1024 * mt - 1 } else gemvCore(F, core) orelse .{ .serial_max = 115200 * mt - 1 },
        // `gbmv.c`, `zgbmv.c`: `m * n < 250000` (125000 complex); their `kl + ku < 15` is the interface's
        .gbmv => .{ .serial_max = if (complex) 124999 else 249999 },
        .ger => if (complex) .{ .serial_max = 36 * s * s * mt } else .{ .serial_max = 2048 * mt, .tiers = comptime &.{upTo(20480 * mt, 3)} },
        // `ztrmv`'s `n*n > 36*s*s*4` and `n*n < 64*s*s*4`, as bounds on `n`
        .trmv => if (complex) .{ .serial_max = 12 * s, .tiers = comptime &.{.{ .below = 16 * s, .max = 2 }} } else .{ .serial_max = 49, .tiers = comptime &.{.{ .below = 500, .max = 2 }} },
        // `syr.c`, `spr.c`, `syr2.c`, `spr2.c`, `zsyr.c`, `zspr.c`: `n < 50`, `n < 500`
        .syr => .{ .serial_max = 49, .tiers = comptime &.{.{ .below = 500, .max = 2 }} },
        // `zher.c`, `zhpr.c`, `zher2.c`, `zhpr2.c`: `n < 100`
        .her => .{ .serial_max = 99 },
        // `symv.c`: `n < 200`; `zsymv.c`, `zhpmv.c`: `n < 50`, `n < 500`
        .symv => if (complex) .{ .serial_max = 49, .tiers = comptime &.{.{ .below = 500, .max = 2 }} } else .{ .serial_max = 199 },
        // `zhemv.c`: `n < MULTI_THREAD_MINIMAL`, 362, then every thread; chunks
        // of `zsymv.c`'s 50 columns, not 361, which ran 2 or 3 workers (#706)
        .hemv => .{ .serial_max = 361, .chunk = 50 },
        // `tbmv.c`, `tpmv.c`, `ztbmv.c`, `ztpmv.c`: `n < 50`, `n < 500`
        .tpmv => .{ .serial_max = 49, .tiers = comptime &.{.{ .below = 500, .max = 2 }} },
        .gemm => (if (complex) null else gemmCore(F, core)) orelse .{ .serial_max = @as(u64, if (complex) 8192 else 65536) * mt, .per_unit = true },
        .syrk => .{ .serial_max = @as(u64, if (complex) 14824 else 109944) * mt },
        // `syr2k.c`: `n * k < 1000`
        .syr2k => .{ .serial_max = 999 },
        .symm => .{ .serial_max = @as(u64, if (complex) 8192 else 65536) * mt },
        .trsm => .{ .serial_max = @as(u64, if (complex) 128 else 256) * mt - 1 },
        // `omatcopy.c`: `rows * cols < 65536 * GEMM_MULTITHREAD_THRESHOLD`
        .omatcopy => .{ .serial_max = 65536 * mt - 1 },
    };
}

/// Whether `core` (`rule`'s) is arm64's `c`.
fn is(core: anytype, comptime c: coretypes.Arm64CoreType) bool {
    const C = @TypeOf(core);
    if (C == coretypes.Arm64CoreType) return core == c;
    if (C == ?coretypes.Arm64CoreType) return if (core) |k| k == c else false;
    return false;
}

/// Whether `core` (`rule`'s) is the x86_64 core `c`.
fn isX86(core: anytype, comptime c: coretypes.X86CoreType) bool {
    const C = @TypeOf(core);
    if (C == coretypes.X86CoreType) return core == c;
    if (C == ?coretypes.X86CoreType) return if (core) |k| k == c else false;
    return false;
}

/// The arm64 cores whose `KERNEL.<core>` takes `?asum_thunderx2t99.c`.
const thunder_asum = [_]coretypes.Arm64CoreType{ .thunderx2t99, .neoversen1, .neoversev1, .neoversev2, .neoversen2, .thunderx3t110, .armv8sve, .a64fx, .armv9sme, .vortex, .vortexm4, .cortexa510, .cortexa710, .cortexx2, .ampere1 };

/// `kernel/arm64/?asum_thunderx2t99.c`: past `n <= 10000`, every thread
/// (issue #728); `?sum_thunderx2t99.c`, on `neoversen1`, `vortex`,
/// `vortexm4` and `ampere1`, alike. Elsewhere `kernel/x86_64/dasum.c`'s.
fn asumCore(core: anytype) ?Rule {
    inline for (thunder_asum) |c| if (is(core, c)) return .{ .serial_max = 10000 };
    return null;
}

/// `kernel/arm64/dot.c`'s tables.
fn dotCore(comptime F: type, core: anytype) ?Rule {
    if (is(core, .neoversev1)) return switch (F) {
        f64 => .{ .serial_max = 64500, .tiers = comptime &.{ upTo(100000, 2), upTo(150000, 4), upTo(260000, 8), upTo(360000, 16), upTo(520000, 24), upTo(1010000, 56) } },
        f32 => .{ .serial_max = 110000, .tiers = comptime &.{ upTo(200000, 2), upTo(280000, 4), upTo(520000, 8), upTo(830000, 16), upTo(1010000, 24) } },
        else => null,
    };
    if (is(core, .a64fx)) return switch (F) {
        f64 => .{ .serial_max = 11000, .tiers = comptime &.{ upTo(20000, 2), upTo(35000, 4), upTo(50000, 6), upTo(440000, 8), upTo(880000, 16), upTo(1020000, 24) } },
        f32 => .{ .serial_max = 22000, .tiers = comptime &.{ upTo(39000, 2), upTo(79000, 4), upTo(120000, 6), upTo(1020000, 8) } },
        else => null,
    };
    return null;
}

/// `interface/gemv.c`'s tables.
fn gemvCore(comptime F: type, core: anytype) ?Rule {
    if (F != f32 and F != f64) return null;
    if (is(core, .neoversev1)) return switch (F) {
        f64 => .{ .serial_max = 8099, .tiers = comptime &.{ .{ .below = 12100, .max = 2 }, .{ .below = 36100, .max = 4 }, .{ .below = 84100, .max = 8 }, .{ .below = 348100, .max = 16 }, .{ .below = 435600, .max = 24 }, .{ .below = 810000, .max = 32 }, .{ .below = 1050625, .max = 40 } } },
        else => .{ .serial_max = 25599, .tiers = comptime &.{ .{ .below = 63001, .max = 4 }, .{ .below = 202500, .max = 8 }, .{ .below = 806404, .max = 16 }, .{ .below = 1638400, .max = 32 } } },
    };
    if (is(core, .neoversev2)) return .{ .serial_max = 24963, .tiers = comptime &.{ .{ .below = 145924, .max = 8 }, .{ .below = 692224, .max = 16 }, .{ .below = 1638400, .max = 32 } } };
    return null;
}

/// `interface/gemm.c`'s tables, single precision only.
fn gemmCore(comptime F: type, core: anytype) ?Rule {
    if (F != f32) return null;
    if (is(core, .neoversev1)) return .{ .serial_max = 262143, .tiers = comptime &.{ .{ .below = 1124864, .max = 6 }, .{ .below = 7880599, .max = 12 }, .{ .below = 17173512, .max = 16 }, .{ .below = 33386248, .max = 20 }, .{ .below = 57066625, .max = 24 }, .{ .below = 91733851, .max = 32 }, .{ .below = 265847707, .max = 40 }, .{ .below = 458314011, .max = 48 }, .{ .below = 729000000, .max = 56 } } };
    if (is(core, .neoversev2)) return .{ .serial_max = 124999, .tiers = comptime &.{ .{ .below = 1092727, .max = 6 }, .{ .below = 2628072, .max = 8 }, .{ .below = 8000000, .max = 12 }, .{ .below = 20346417, .max = 16 }, .{ .below = 57066625, .max = 24 }, .{ .below = 91125000, .max = 28 }, .{ .below = 238328000, .max = 40 }, .{ .below = 454756609, .max = 48 }, .{ .below = 857375000, .max = 56 }, .{ .below = 1073741824, .max = 64 } } };
    return null;
}

/// `dims` multiplied, saturating: a rule's `work`.
pub fn product(dims: []const u64) u64 {
    var p: u64 = 1;
    for (dims) |d| p *|= d;
    return p;
}
