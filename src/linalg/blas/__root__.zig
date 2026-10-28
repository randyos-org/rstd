//! Derived from code in OpenBLAS

const std = @import("std");
const sysinfo = @import("builtin");
const portable = @import("portable.zig");
const threading = @import("threading.zig");
const apple_accel = @import("apple_accel.zig");

pub const types = @import("types.zig");
pub const cpu = @import("cpu/__root__.zig");
pub const Vtable = @import("Vtable.zig");
pub const Config = @import("Config.zig");
/// The `@Vector` kernels some cores install; `tools/blasbench` races them.
pub const vector = @import("vector.zig");
const blas1 = @import("blas1/__root__.zig");
const blas2 = @import("blas2/__root__.zig");
const blas3 = @import("blas3/__root__.zig");
const ext = @import("ext/__root__.zig");
const l1 = blas1.interface;
const l2 = blas2.interface;
const tri = blas2.triangular;
const band = blas2.banded;
const l3 = blas3.interface;
// pub const util = @import("util.zig");
const util = @import("util.zig");

// ---- the BLAS routines ---------------------------------------------------------
//
// Fortran argument order, column-major, sizes and strides `isize`, the element
// type explicit. Each takes `io` after its allocator, or first without one
// (`threading.md`); bad arguments return named errors in `xerbla`'s order.

pub const Trans = types.Trans;
pub const Uplo = types.Uplo;
pub const Diag = types.Diag;
pub const Side = types.Side;
pub const Extreme = l1.Extreme;
pub const HalfGemv = l2.HalfGemv;

/// Level 1: a stride past the address space, or `error.Canceled`.
pub const Level1Error = l1.Error;
/// Level 2: bad arguments, `error.OutOfMemory` or `error.Canceled`.
pub const Level2Error = l2.Error;
/// `gemvHalf`, which allocates nothing: `Level2Error` less `OutOfMemory`.
pub const GemvHalfError = l2.HalfError;
/// Level 3: bad arguments, `error.OutOfMemory` or `error.Canceled`.
pub const Level3Error = l3.Error;

pub const copy = l1.copy;
pub const swap = l1.swap;
pub const scal = l1.scal;
pub const scalReal = l1.scalReal;
pub const dot = l1.dot;
pub const dotc = l1.dotc;
pub const dsdot = l1.dsdot;
pub const sdsdot = l1.sdsdot;
pub const nrm2 = l1.nrm2;
pub const asum = l1.asum;
pub const sum = l1.sum;
pub const index = l1.index;
pub const extreme = l1.extreme;
pub const rot = l1.rot;
pub const rotm = l1.rotm;
pub const axpy = l1.axpy;
pub const axpyc = l1.axpyc;
pub const axpby = l1.axpby;
pub const sbdot = l1.sbdot;
pub const bscal = l1.bscal;
pub const sbtobf16 = l1.sbtobf16;
pub const bf16to = l1.bf16to;
pub const rotg = l1.rotg;
pub const zrotg = l1.zrotg;
pub const rotmg = l1.rotmg;
pub const cabs1 = l1.cabs1;

pub const gemv = l2.gemv;
pub const gemvHalf = l2.gemvHalf;
pub const ger = l2.ger;
pub const gerc = l2.gerc;
pub const symv = l2.symv;
pub const hemv = l2.hemv;
pub const trmv = tri.trmv;
pub const trsv = tri.trsv;
pub const tbmv = tri.tbmv;
pub const tbsv = tri.tbsv;
pub const tpmv = tri.tpmv;
pub const tpsv = tri.tpsv;
pub const gbmv = band.gbmv;
pub const sbmv = band.sbmv;
pub const hbmv = band.hbmv;
pub const spmv = band.spmv;
pub const hpmv = band.hpmv;
pub const syr = band.syr;
pub const spr = band.spr;
pub const her = band.her;
pub const hpr = band.hpr;
pub const syr2 = band.syr2;
pub const spr2 = band.spr2;
pub const her2 = band.her2;
pub const hpr2 = band.hpr2;

pub const gemm = l3.gemm;
pub const symm = l3.symm;
pub const hemm = l3.hemm;
pub const syrk = l3.syrk;
pub const herk = l3.herk;
pub const syr2k = l3.syr2k;
pub const her2k = l3.her2k;
pub const trmm = l3.trmm;
pub const trsm = l3.trsm;

// ---- selection and the vtable ----------------------------------------------------

const InitPhase = enum(u8) { uninit, running, ready };

pub var config: Config = .{};
/// Null until `init` selects a core (and under a test's own table).
var core_type: ?cpu.CoreType = null;
var _vtable: Vtable = generic_vtable;

/// The vtable before `init` and after `shutdown`: the portable kernels, where
/// there are any.
const generic_vtable: Vtable = blk: {
    @setEvalBranchQuota(100_000);
    var vt: Vtable = .{};
    portable.install(&vt);
    break :blk vt;
};

/// A vtable slot's name.
pub const Slot = util.Slot;

var phase = std.atomic.Value(InitPhase).init(.uninit);
/// How many times a core was selected; the tests watch the handshake with it.
var inits = std.atomic.Value(u32).init(0);

/// The vtable, initialized first if nothing has done so (or after `shutdown`).
/// A routine fetches it once per call and passes it down, or uses `kernel`
/// once; code under such a call reads `table` instead.
pub inline fn vtable() *const Vtable {
    if (phase.load(.acquire) != .ready) initCold();
    return &_vtable;
}

fn initCold() void {
    @branchHint(.cold);
    _ = initWith(processEnv());
}

/// The vtable, for code a routine reaches after its own `vtable()` call.
pub inline fn table() *const Vtable {
    if (comptime sysinfo.optimize.runtimeSafety()) std.debug.assert(phase.load(.monotonic) == .ready);
    return &_vtable;
}

/// The kernel in vtable slot `slot` (`.daxpy_k`), as the reference's `DAXPY_K`
/// is `gotoblas->daxpy_k`. A slot with no kernel for the running core, and no
/// portable one, is a bug in the caller.
pub inline fn kernel(comptime slot: Slot) @typeInfo(@FieldType(Vtable, @tagName(slot))).optional.child {
    return @field(vtable(), @tagName(slot)).?;
}

/// Selects the core and fills the vtable with its kernels, once; true if this
/// call did it. The first kernel call does the same on its own, reading the
/// process environment where a library can (`processEnv`); this call is how
/// `OPENBLAS_*` is honored elsewhere, and it must come before the first kernel
/// call to count. Later calls do nothing until `shutdown`. The values kept are
/// copied, so `env` may go away after.
pub fn init(env: *const std.process.Environ.Map) bool {
    return initWith(.{ .map = env });
}

/// A mutex would drag in threading support we do not want on freestanding
/// targets, so this is a three-state handshake: the winner of the
/// compare-exchange initializes while any other thread spins until it publishes.
fn initWith(env: Config.Env) bool {
    if (!claim()) return false;
    config = .fromEnv(env);
    config.cpu_count = threading.cpuCount();

    const chosen = cpu.choose(config.coreType(), config.ld_hwcap_mask, LogSink{});
    // the core whose kernels run, as `gotoblas_corename` names it
    const ct = chosen.served orelse chosen.core;
    core_type = ct;
    blasWarn(2, "core: {s}", .{@tagName(ct)});
    if (comptime sysinfo.target.cpu.arch.isPowerpc()) {
        if (config.get_openblas_coretype and comptime sysinfo.target.os.tag != .freestanding)
            std.debug.print("Core: {s}\n", .{@tagName(ct)});
    }

    _vtable = cpu.vtableFromCoreType(ct);
    if (_vtable.init) |f| f();
    apple_accel.start(config.accelerate);
    _ = inits.fetchAdd(1, .monotonic);
    phase.store(.ready, .release);
    return true;
}

/// Moves `uninit` to `running` for this caller; false once another caller has
/// published the table.
fn claim() bool {
    while (true) {
        switch (phase.load(.acquire)) {
            .ready => return false,
            .running => std.atomic.spinLoopHint(),
            .uninit => if (phase.cmpxchgWeak(.uninit, .running, .acquire, .monotonic) == null) return true,
        }
    }
}

/// The process environment where a library can read it: the Windows PEB, or
/// libc's `environ` (`_NSGetEnviron()` on Darwin, where a dylib has no
/// `environ`). Zig hands the block to `main` only, so elsewhere there is none.
fn processEnv() Config.Env {
    if (comptime sysinfo.target.os.tag == .windows) return .windows;
    if (comptime !sysinfo.link_libc) return .none;
    const envp: [*:null]const ?[*:0]const u8 = switch (sysinfo.target.os.tag) {
        .freestanding, .other, .uefi => return .none,
        else => if (comptime sysinfo.target.os.tag.isDarwin()) _NSGetEnviron().* else std.c.environ,
    };
    return .{ .posix = std.mem.span(envp) };
}

extern "c" fn _NSGetEnviron() *[*:null]?[*:0]u8;

/// `cpu.choose`'s log: `OPENBLAS_VERBOSE` decides what shows.
const LogSink = struct {
    pub fn warn(_: LogSink, level: i32, comptime fmt: []const u8, args: anytype) void {
        blasWarn(level, fmt, args);
    }
};

/// A diagnostic, shown when `OPENBLAS_VERBOSE` is at least `level`. Only where
/// there is somewhere to print to.
fn blasWarn(level: i32, comptime fmt: []const u8, args: anytype) void {
    if (comptime sysinfo.target.os.tag == .freestanding) return;
    if (config.verbose < level) return;
    std.debug.print("OpenBLAS: " ++ fmt ++ "\n", args);
}

/// Back to the portable kernels and the unset state; the next kernel call, or
/// `init`, selects a core again. No kernel call may be running (the
/// reference's `gotoblas_quit` has the same rule).
pub fn shutdown() void {
    if (phase.cmpxchgStrong(.ready, .running, .acquire, .monotonic) != null) return;
    // if (vtable.quit) |f| f();
    _ = apple_accel.route(false);
    _vtable = generic_vtable;
    config = .{};
    core_type = null;
    phase.store(.uninit, .release);
}

/// Turns the routines' use of Apple Accelerate on or off (macOS, `s`, `d`, `c`
/// and `z`), after initializing if need be; true if they use it now. It is
/// loaded at the first init where `RSTD_BLAS_ACCELERATE` is not 0 and never
/// unloaded; before that, this cannot turn it on.
pub fn useAccelerate(enable: bool) bool {
    _ = vtable();
    return apple_accel.route(enable);
}

/// The selected core, selecting it first if need be; null under a test's own
/// table (`testUse`).
pub fn coreType() ?cpu.CoreType {
    _ = vtable();
    return core_type;
}

pub fn coreName() []const u8 {
    return if (coreType()) |ct| @tagName(ct) else "test";
}

/// What `testUse` replaced.
pub const TestSaved = struct { vt: Vtable, ready: bool, accelerate: bool };

/// Puts `vt` in place of the selected kernels, without selecting a core
/// first, until `testRestore`; the routines run them, not Accelerate.
pub fn testUse(vt: *const Vtable) TestSaved {
    if (!sysinfo.is_test) @compileError("tests only");
    const saved: TestSaved = .{ .vt = _vtable, .ready = !claim(), .accelerate = apple_accel.routing() };
    _ = apple_accel.route(false);
    _vtable = vt.*;
    phase.store(.ready, .release);
    return saved;
}

/// The selected table itself, for a test that adjusts one parameter in place
/// (after `testUse`, so a core is not selected under it).
pub fn testTable() *Vtable {
    if (!sysinfo.is_test) @compileError("tests only");
    return &_vtable;
}

pub fn testRestore(saved: TestSaved) void {
    if (!sysinfo.is_test) @compileError("tests only");
    _vtable = saved.vt;
    _ = apple_accel.route(saved.accelerate);
    if (!saved.ready) phase.store(.uninit, .release);
}

const zigstd = @import("../../zigstd.zig");

// lives here, not in `Vtable.zig`, so that file stays import-free for the
// per-core kernel objects, which use its type (see `kernels/__root__.zig`)
test "VTable default-initializes with every kernel pointer null" {
    const vt: Vtable = .{};
    try zigstd.testExpectTrue(vt.sgemm_kernel == null);
    try zigstd.testExpectTrue(vt.dgemm_kernel == null);
    try zigstd.testExpectTrue(vt.cgemm_kernel_n == null);
    try zigstd.testExpectTrue(vt.zgemm_kernel_n == null);
    try zigstd.testExpectTrue(vt.init == null);
    try zigstd.testExpectEqual(@as(i32, 0), vt.dtb_entries);
}

test "a kernel call initializes by itself, and init after that does nothing" {
    shutdown();
    defer _ = vtable();
    const before = inits.load(.monotonic);
    try zigstd.testExpectTrue(phase.load(.acquire) == .uninit);
    try zigstd.testExpectTrue(kernel(.daxpy_k) == table().daxpy_k.?);
    try zigstd.testExpectTrue(phase.load(.acquire) == .ready);
    try zigstd.testExpectEqual(before + 1, inits.load(.monotonic));
    try zigstd.testExpectTrue(table().dtb_entries > 0);
    try zigstd.testExpectTrue(coreType() != null);
    var env: std.process.Environ.Map = .init(zigstd.testAllocator);
    defer env.deinit();
    try zigstd.testExpectTrue(!init(&env));
    try zigstd.testExpectEqual(before + 1, inits.load(.monotonic));
}

test "init keeps its own copy of the environment, shutdown goes back to the portable kernels" {
    shutdown();
    defer _ = vtable();
    try zigstd.testExpectTrue(_vtable.daxpy_k == generic_vtable.daxpy_k);
    const detected = @tagName(cpu.detect());
    // the core that runs it: a detected core not linked into this build falls back
    const chosen = cpu.choose(detected, null, LogSink{});
    const runs = @tagName(chosen.served orelse chosen.core);
    {
        var env: std.process.Environ.Map = .init(zigstd.testAllocator);
        defer env.deinit();
        try env.put("OPENBLAS_CORETYPE", detected);
        try env.put("OPENBLAS_NUM_THREADS", "3");
        try zigstd.testExpectTrue(init(&env));
    }
    // the map is gone; the copy is not
    try zigstd.testExpectEqualStrings(detected, config.coreType().?);
    try zigstd.testExpectEqual(@as(i32, 3), config.num_threads);
    try zigstd.testExpectEqualStrings(runs, coreName());

    shutdown();
    try zigstd.testExpectEqual(@as(i32, 0), config.num_threads);
    try zigstd.testExpectTrue(_vtable.daxpy_k == generic_vtable.daxpy_k);
    var env: std.process.Environ.Map = .init(zigstd.testAllocator);
    defer env.deinit();
    try env.put("OPENBLAS_CORETYPE", "no-such-core");
    try zigstd.testExpectTrue(init(&env));
    try zigstd.testExpectEqualStrings(runs, coreName());
}

test "a test's table needs no core, and goes back to the unset state" {
    shutdown();
    defer _ = vtable();
    const before = inits.load(.monotonic);
    const vt: Vtable = generic_vtable;
    const saved = testUse(&vt);
    try zigstd.testExpectTrue(phase.load(.acquire) == .ready);
    try zigstd.testExpectTrue(coreType() == null);
    testRestore(saved);
    try zigstd.testExpectTrue(phase.load(.acquire) == .uninit);
    try zigstd.testExpectEqual(before, inits.load(.monotonic));
}

var race_seen: [4]?cpu.CoreType = @splat(null);

fn raceKernel(i: usize) void {
    race_seen[i] = if (kernel(.daxpy_k) == _vtable.daxpy_k.?) core_type else null;
}

test "threads racing the first kernel call select the core once and all see it" {
    if (sysinfo.single_threaded) return error.SkipZigTest;
    shutdown();
    defer _ = vtable();
    race_seen = @splat(null);
    const before = inits.load(.monotonic);
    var threads: [race_seen.len]std.Thread = undefined;
    {
        // a failed spawn still joins the threads already running
        var spawned: usize = 0;
        defer for (threads[0..spawned]) |t| t.join();
        for (&threads, 0..) |*t, i| {
            t.* = try std.Thread.spawn(.{}, raceKernel, .{i});
            spawned += 1;
        }
    }
    try zigstd.testExpectEqual(before + 1, inits.load(.monotonic));
    try zigstd.testExpectTrue(phase.load(.acquire) == .ready);
    for (race_seen) |seen| try zigstd.testExpectTrue(seen != null and seen.? == core_type.?);
    try zigstd.testExpectTrue(_vtable.dtb_entries > 0);
}

/// The `q`/`x` kernel slots' names.
const qx_slots = blk: {
    @setEvalBranchQuota(100_000);
    const info = @typeInfo(Vtable).@"struct";
    var names: []const []const u8 = &.{};
    for (info.field_names, info.field_types) |f, F| {
        const qx = f[0] == 'q' or f[0] == 'x' or (f[0] == 'i' and (f[1] == 'q' or f[1] == 'x'));
        if (qx and @typeInfo(F) == .optional) names = names ++ .{f};
    }
    break :blk names;
};

test "no core puts a kernel of its own in a q or x slot" {
    @setEvalBranchQuota(100_000);
    const registry = @import("asm_registry.zig");
    var base: Vtable = .{};
    portable.install(&base);
    inline for (@typeInfo(cpu.CoreType).@"enum".field_names) |name| {
        if (registry.isLinked(@field(cpu.CoreType, name))) try qxUntouched(&base, @field(cpu.CoreType, name));
    }
}

/// One core per call: inlined into the loop, each core's table and install
/// would keep a frame slot of their own (~400 KB on a 64-bit debug build).
fn qxUntouched(base: *const Vtable, comptime core: cpu.CoreType) !void {
    var vt = base.*;
    cpu.kernel.install(&vt, core);
    inline for (qx_slots) |f| try zigstd.testExpectTrue(@field(vt, f) == @field(base, f));
}

test {
    _ = @import("abi.zig");
    _ = @import("apple_accel.zig");
    _ = @import("threading.zig");
    _ = @import("_test/__root__.zig");
    zigstd.testRefAllDecls(@This());
}
