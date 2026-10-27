//! Everything the reference reads from the environment at init, in one
//! place. Defaults match the reference's "unset" behavior (0 / absent).

const Self = @This();

const std = @import("std");
const sysinfo = @import("builtin");
const zigstd = @import("../../zigstd.zig");

/// `OPENBLAS_VERBOSE` -- 0 quiet, 2 prints the selected core name.
verbose: i32 = 0,
/// `OPENBLAS_BLOCK_FACTOR` -- scales the GEMM blocking parameters.
block_factor: i32 = 0,
/// `OPENBLAS_THREAD_TIMEOUT` -- thread-pool spin timeout, log2 ticks.
thread_timeout: u32 = 0,
/// `OPENBLAS_NUM_THREADS`, falling back to
/// `OPENBLAS_DEFAULT_NUM_THREADS`. The reference lets the former win
/// unless it is unset and the latter is not.
num_threads: i32 = 0,
/// `GOTO_NUM_THREADS` -- legacy GotoBLAS spelling.
goto_num_threads: i32 = 0,
/// `OMP_NUM_THREADS` -- consulted only in OpenMP builds.
omp_num_threads: i32 = 0,
/// `OMP_ADAPTIVE`.
omp_adaptive: i32 = 0,
/// `OPENBLAS_CORETYPE` -- forces a core type by name, bypassing CPU
/// detection entirely. Invaluable for benchmarking and for working
/// around a bad detection result without a rebuild. Read in every
/// `dynamic_*.c`, not just x86's. Kept here, not borrowed from the
/// environment, which may change or go away; see `coreType`.
core_type_buf: [value_max]u8 = @splat(0),
core_type_len: u8 = 0,

/// `GET_OPENBLAS_CORETYPE` -- POWER-specific (`dynamic_power.c` only).
/// A bare presence check (`if (getenv(...))`), not a parsed value: if
/// set to anything at all, including an empty string, the selected
/// core name is echoed straight to stderr in addition to the normal
/// `openblas_warning` mechanism.
get_openblas_coretype: bool = false,

/// `GOTO_DEBUG` -- legacy per-call debug print, gating the
/// `PRINT_DEBUG_CNAME`/`PRINT_DEBUG_NAME` macros (`common.h`).
/// Does *not* clamp negative values to 0.
goto_debug: i32 = 0,

/// `OPENBLAS_MAIN_FREE` -- whether the main thread's NUMA
/// shared-memory buffer gets freed (`driver/others/init.c`).
/// Checks `"OPENBLAS_MAIN_FREE" || "GOTOBLAS_MAIN_FREE"`, collapsed to one bool.
main_free: bool = false,

/// `OPENBLAS_L2_SIZE` -- overrides the detected L2 cache size used by
/// `kernel/setparam-ref.c`'s tuning heuristics.
/// No negative clamp); 0 means "not overridden"
/// (an L2 size of 0 is never a real value to tune for).
l2_size_override: i32 = 0,

/// Not from the environment: the CPU count `blas.init` caches
/// (`threading.cpuCount`), 1 until then.
cpu_count: u32 = 1,

/// `LD_HWCAP_MASK` -- s390x-only (`cpuid_zarch.h`). ANDed against
/// `getauxval(AT_HWCAP)` to let the user mask *out* specific feature
/// bits the kernel reports, even if the hardware genuinely has them.
/// Parsed with `strtoul(p, NULL, 0)` -- base-0, so `0x`/`0`-prefixed
/// values pick hex/octal the way the C source does.
ld_hwcap_mask: ?u64 = null,

/// `RSTD_BLAS_ACCELERATE` (rstd's own) -- 0 keeps the routines off Apple
/// Accelerate on macOS, which is then never loaded; anything else, or unset,
/// uses it (`apple_accel.zig`).
accelerate: bool = true,

/// The longest value read; longer ones read as unparseable (a NUL), so they
/// fall back to the default (or, for a core name, to detection).
pub const value_max = 64;

/// Where the variables come from: a map the program passes (`blas.init`), or
/// the process's own block, read one variable at a time with no allocation.
pub const Env = union(enum) {
    map: *const std.process.Environ.Map,
    /// `NAME=value` strings, as `environ` holds them.
    posix: []const ?[*:0]const u8,
    /// The process's, through the PEB.
    windows,
    none,

    /// `key`'s value, empty if set to nothing; copied into `buf` where it
    /// cannot be borrowed (Windows).
    pub fn get(e: Env, comptime key: []const u8, buf: *[value_max]u8) ?[]const u8 {
        const v: []const u8 = switch (e) {
            .map => |m| m.get(key) orelse return null,
            .posix => |entries| for (entries) |opt| {
                const entry = std.mem.span(opt orelse continue);
                if (entry.len > key.len and entry[key.len] == '=' and std.mem.eql(u8, entry[0..key.len], key)) break entry[key.len + 1 ..];
            } else return null,
            .windows => if (comptime sysinfo.target.os.tag == .windows) blk: {
                const w = std.process.Environ.getWindows(.{ .block = .global }, std.unicode.wtf8ToWtf16LeStringLiteral(key)) orelse return null;
                if (w.len * 3 > buf.len) break :blk tooLong(buf);
                break :blk buf[0..std.unicode.wtf16LeToWtf8(buf, w)];
            } else return null,
            .none => return null,
        };
        return if (v.len > value_max) tooLong(buf) else v;
    }

    fn tooLong(buf: *[value_max]u8) []const u8 {
        buf[0] = 0;
        return buf[0..1];
    }
};

fn asEnv(env: anytype) Env {
    return if (@TypeOf(env) == Env) env else .{ .map = env };
}

pub fn init(env: *const std.process.Environ.Map) Self {
    return fromEnv(.{ .map = env });
}

pub fn fromEnv(env: Env) Self {
    var c: Self = .{};
    c.verbose = readInt(env, "OPENBLAS_VERBOSE");
    c.block_factor = readInt(env, "OPENBLAS_BLOCK_FACTOR");
    c.thread_timeout = @intCast(readInt(env, "OPENBLAS_THREAD_TIMEOUT"));

    const default_threads = readInt(env, "OPENBLAS_DEFAULT_NUM_THREADS");
    const threads = readInt(env, "OPENBLAS_NUM_THREADS");
    c.num_threads = if (threads != 0 or default_threads == 0) threads else default_threads;

    c.goto_num_threads = readInt(env, "GOTO_NUM_THREADS");
    c.omp_num_threads = readInt(env, "OMP_NUM_THREADS");
    c.omp_adaptive = readInt(env, "OMP_ADAPTIVE");
    var buf: [value_max]u8 = undefined;
    if (readValue(env, "OPENBLAS_CORETYPE", &buf)) |name| {
        @memcpy(c.core_type_buf[0..name.len], name);
        c.core_type_len = @intCast(name.len);
    }

    c.get_openblas_coretype = isSet(env, "GET_OPENBLAS_CORETYPE");
    c.goto_debug = readIntRaw(i32, env, "GOTO_DEBUG", 0) orelse 0;
    c.main_free = (readIntRaw(i32, env, "OPENBLAS_MAIN_FREE", 0) orelse 0) != 0 or
        (readIntRaw(i32, env, "GOTOBLAS_MAIN_FREE", 0) orelse 0) != 0;
    c.l2_size_override = readIntRaw(i32, env, "OPENBLAS_L2_SIZE", 0) orelse 0;
    c.ld_hwcap_mask = readIntRaw(u64, env, "LD_HWCAP_MASK", null);
    c.accelerate = (readIntRaw(i32, env, "RSTD_BLAS_ACCELERATE", 1) orelse 1) != 0;
    return c;
}

/// `OPENBLAS_CORETYPE`, if set.
pub fn coreType(c: *const Self) ?[]const u8 {
    return if (c.core_type_len == 0) null else c.core_type_buf[0..c.core_type_len];
}

/// The thread ceiling: `num_threads`, else `goto_num_threads`; null for none.
pub fn threadCeiling(c: *const Self) ?u32 {
    if (c.num_threads > 0) return @intCast(c.num_threads);
    if (c.goto_num_threads > 0) return @intCast(c.goto_num_threads);
    return null;
}

/// Reads `key` into `buf` (or borrows it), returning null when unset or empty.
pub fn readValue(env: anytype, comptime key: []const u8, buf: *[value_max]u8) ?[]const u8 {
    const v = asEnv(env).get(key, buf) orelse return null;
    return if (v.len == 0) null else v;
}

/// clamps negatives to 0.
pub fn readInt(env: anytype, comptime key: []const u8) i32 {
    const parsed = readIntRaw(i32, env, key, 0) orelse 0;
    return if (parsed < 0) 0 else parsed;
}

/// clamps unset to `default`.
pub fn readIntRaw(comptime T: type, env: anytype, comptime key: []const u8, default: ?T) ?T {
    var buf: [value_max]u8 = undefined;
    const v = readValue(env, key, &buf) orelse return default;
    return std.fmt.parseInt(T, std.mem.trim(u8, v, " \t\r\n"), 0) catch default;
}

/// True if `key` is set in the environment at all, regardless of value --
/// even an empty string counts, matching a bare C `if (getenv(...))`.
/// `readValue` treats an empty value as "not set," so this can't just
/// delegate to it.
pub fn isSet(env: anytype, comptime key: []const u8) bool {
    var buf: [value_max]u8 = undefined;
    return asEnv(env).get(key, &buf) != null;
}

/// Builds an empty `Environ.Map` for a test, pre-populated with `pairs`.
fn testEnv(pairs: []const [2][]const u8) std.process.Environ.Map {
    var env: std.process.Environ.Map = .init(zigstd.testAllocator);
    for (pairs) |kv| env.put(kv[0], kv[1]) catch unreachable;
    return env;
}

test "unset variables yield every default" {
    var env = testEnv(&.{});
    defer env.deinit();
    var buf: [value_max]u8 = undefined;

    try zigstd.testExpectEqual(@as(i32, 0), readInt(&env, "OPENBLAS_DEFINITELY_NOT_SET_XYZ"));
    try zigstd.testExpectTrue(readValue(&env, "OPENBLAS_DEFINITELY_NOT_SET_XYZ", &buf) == null);
    try zigstd.testExpectEqual(@as(?i32, 0), readIntRaw(i32, &env, "OPENBLAS_DEFINITELY_NOT_SET_XYZ", 0));
    try zigstd.testExpectTrue(!isSet(&env, "OPENBLAS_DEFINITELY_NOT_SET_XYZ"));

    const c = Self.init(&env);
    try zigstd.testExpectEqual(Self{}, c);
}

test "readInt clamps a negative value to 0" {
    var env = testEnv(&.{.{ "OPENBLAS_BLOCK_FACTOR", "-5" }});
    defer env.deinit();

    try zigstd.testExpectEqual(@as(i32, 0), readInt(&env, "OPENBLAS_BLOCK_FACTOR"));
}

test "readIntRaw does not clamp negatives, unlike readInt" {
    var env = testEnv(&.{.{ "GOTO_DEBUG", "-5" }});
    defer env.deinit();

    try zigstd.testExpectEqual(@as(?i32, -5), readIntRaw(i32, &env, "GOTO_DEBUG", 0));
    try zigstd.testExpectEqual(@as(i32, 0), readInt(&env, "GOTO_DEBUG"));

    const c = Self.init(&env);
    try zigstd.testExpectEqual(@as(i32, -5), c.goto_debug);
}

test "a malformed value falls back to the default instead of erroring" {
    var env = testEnv(&.{.{ "OPENBLAS_VERBOSE", "not-a-number" }});
    defer env.deinit();

    try zigstd.testExpectEqual(@as(?i32, 0), readIntRaw(i32, &env, "OPENBLAS_VERBOSE", 0));
}

test "OPENBLAS_NUM_THREADS wins over OPENBLAS_DEFAULT_NUM_THREADS when both are set" {
    var env = testEnv(&.{
        .{ "OPENBLAS_NUM_THREADS", "4" },
        .{ "OPENBLAS_DEFAULT_NUM_THREADS", "8" },
    });
    defer env.deinit();

    try zigstd.testExpectEqual(@as(i32, 4), Self.init(&env).num_threads);
}

test "OPENBLAS_DEFAULT_NUM_THREADS is used only when OPENBLAS_NUM_THREADS is unset" {
    var env = testEnv(&.{.{ "OPENBLAS_DEFAULT_NUM_THREADS", "8" }});
    defer env.deinit();

    try zigstd.testExpectEqual(@as(i32, 8), Self.init(&env).num_threads);
}

test "neither thread-count variable set yields 0 threads" {
    var env = testEnv(&.{});
    defer env.deinit();

    try zigstd.testExpectEqual(@as(i32, 0), Self.init(&env).num_threads);
}

test "the thread ceiling is OPENBLAS_NUM_THREADS, else GOTO_NUM_THREADS" {
    var env = testEnv(&.{.{ "GOTO_NUM_THREADS", "6" }});
    defer env.deinit();
    try zigstd.testExpectEqual(@as(?u32, 6), Self.init(&env).threadCeiling());
    try env.put("OPENBLAS_NUM_THREADS", "3");
    try zigstd.testExpectEqual(@as(?u32, 3), Self.init(&env).threadCeiling());
    try env.put("OPENBLAS_NUM_THREADS", "-2");
    try env.put("GOTO_NUM_THREADS", "0");
    try zigstd.testExpectEqual(@as(?u32, null), Self.init(&env).threadCeiling());
}

test "GET_OPENBLAS_CORETYPE is a bare-presence check, even set to empty" {
    var env = testEnv(&.{.{ "GET_OPENBLAS_CORETYPE", "" }});
    defer env.deinit();

    try zigstd.testExpectTrue(isSet(&env, "GET_OPENBLAS_CORETYPE"));
    var buf: [value_max]u8 = undefined;
    try zigstd.testExpectTrue(readValue(&env, "GET_OPENBLAS_CORETYPE", &buf) == null);
    try zigstd.testExpectTrue(Self.init(&env).get_openblas_coretype);
}

test "OPENBLAS_MAIN_FREE and GOTOBLAS_MAIN_FREE are either-or" {
    {
        var env = testEnv(&.{.{ "OPENBLAS_MAIN_FREE", "1" }});
        defer env.deinit();
        try zigstd.testExpectTrue(Self.init(&env).main_free);
    }
    {
        var env = testEnv(&.{.{ "GOTOBLAS_MAIN_FREE", "1" }});
        defer env.deinit();
        try zigstd.testExpectTrue(Self.init(&env).main_free);
    }
    {
        var env = testEnv(&.{});
        defer env.deinit();
        try zigstd.testExpectTrue(!Self.init(&env).main_free);
    }
}

test "LD_HWCAP_MASK base-0 parsing: matches strtoul for hex, diverges for bare-leading-zero octal" {
    {
        var env = testEnv(&.{.{ "LD_HWCAP_MASK", "0xff" }});
        defer env.deinit();
        try zigstd.testExpectEqual(@as(?u64, 0xff), readIntRaw(u64, &env, "LD_HWCAP_MASK", null));
    }
    {
        // C's strtoul(p, NULL, 0) would read "010" as octal 8; Zig's
        // parseInt with base 0 does not infer octal from a bare leading
        // zero, only from an explicit "0o" prefix, so this parses as
        // decimal 10, not 8.
        var env = testEnv(&.{.{ "LD_HWCAP_MASK", "010" }});
        defer env.deinit();
        try zigstd.testExpectEqual(@as(?u64, 10), readIntRaw(u64, &env, "LD_HWCAP_MASK", null));
    }
    {
        var env = testEnv(&.{.{ "LD_HWCAP_MASK", "0o10" }});
        defer env.deinit();
        try zigstd.testExpectEqual(@as(?u64, 8), readIntRaw(u64, &env, "LD_HWCAP_MASK", null));
    }
}

test "RSTD_BLAS_ACCELERATE=0 turns Accelerate off; unset or anything else leaves it on" {
    for ([_]struct { v: ?[]const u8, on: bool }{
        .{ .v = null, .on = true },
        .{ .v = "0", .on = false },
        .{ .v = " 0 ", .on = false },
        .{ .v = "1", .on = true },
        .{ .v = "no", .on = true },
    }) |c| {
        var env = testEnv(&.{});
        defer env.deinit();
        if (c.v) |v| try env.put("RSTD_BLAS_ACCELERATE", v);
        try zigstd.testExpectEqual(c.on, Self.init(&env).accelerate);
    }
}

test "a fully populated environment produces a matching config" {
    var env = testEnv(&.{
        .{ "OPENBLAS_VERBOSE", "2" },
        .{ "OPENBLAS_CORETYPE", "Haswell" },
        .{ "OPENBLAS_L2_SIZE", "262144" },
    });
    defer env.deinit();

    const c = Self.init(&env);
    env.deinit();
    env = testEnv(&.{});
    try zigstd.testExpectEqual(@as(i32, 2), c.verbose);
    try zigstd.testExpectEqualStrings("Haswell", c.coreType().?);
    try zigstd.testExpectEqual(@as(i32, 262144), c.l2_size_override);
}

test "the process block is read one variable at a time, and over-long values read as unparseable" {
    const xs: [80:0]u8 = @splat('x');
    const ones: [70:0]u8 = @splat('1');
    const entries = [_]?[*:0]const u8{ "OPENBLAS_VERBOSEX=9", "OPENBLAS_VERBOSE=3", "OPENBLAS_CORETYPE=" ++ xs, "OPENBLAS_NUM_THREADS=" ++ ones, null };
    const c = fromEnv(.{ .posix = entries[0..4] });
    try zigstd.testExpectEqual(@as(i32, 3), c.verbose);
    try zigstd.testExpectEqual(@as(i32, 0), c.num_threads);
    try zigstd.testExpectEqualStrings("\x00", c.coreType().?);
}
