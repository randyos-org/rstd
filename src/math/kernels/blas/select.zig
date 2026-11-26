//! Core selection and the vtable: the first kernel call, or `init`, picks
//! the core and installs its kernels, once; `shutdown` goes back. The
//! `__root__.zig` API re-exports it; the state is `pub` for `_test/select.zig`.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("__root__.zig");
const portable = @import("portable.zig");
const threading = @import("threading.zig");
const apple_accel = @import("apple_accel.zig");
const cpu = @import("cpu/__root__.zig");
const util = @import("util.zig");
const Config = @import("Config.zig");
const Vtable = @import("Vtable.zig");

const InitPhase = enum(u8) { uninit, running, ready };

/// Null until `init` selects a core (and under a test's own table).
pub var core_type: ?cpu.CoreType = null;
pub var _vtable: Vtable = generic_vtable;

/// The vtable before `init` and after `shutdown`: the portable kernels, where
/// there are any.
pub const generic_vtable: Vtable = blk: {
    @setEvalBranchQuota(1_000_000);
    var vt: Vtable = .{};
    portable.install(&vt);
    break :blk vt;
};

/// A vtable slot's name.
pub const Slot = util.Slot;

pub var phase = std.atomic.Value(InitPhase).init(.uninit);
/// How many times a core was selected; the tests watch the handshake with it.
pub var inits = std.atomic.Value(u32).init(0);

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
/// Asserts, in safe builds, that it is installed: reaching a kernel first is
/// a broken invariant.
pub inline fn table() *const Vtable {
    if (comptime sysinfo.optimize.runtimeSafety()) std.debug.assert(phase.load(.monotonic) == .ready);
    return &_vtable;
}

/// The vtable for `T`: a type without slots runs the portable kernels and
/// never reads it.
pub fn tableOf(comptime T: type) *const Vtable {
    return if (comptime portable.letter(T) != null) vtable() else &no_table;
}

const no_table: Vtable = .{};

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
    blas.config = .fromEnv(env);
    blas.config.cpu_count = threading.cpuCount();

    const chosen = cpu.choose(blas.config.coreType(), blas.config.ld_hwcap_mask, LogSink{});
    // the core whose kernels run, as `gotoblas_corename` names it
    const ct = chosen.served orelse chosen.core;
    core_type = ct;
    blasWarn(2, "core: {s}", .{@tagName(ct)});
    if (comptime sysinfo.target.cpu.arch.isPowerpc()) {
        if (blas.config.get_openblas_coretype and comptime sysinfo.target.os.tag != .freestanding)
            std.debug.print("Core: {s}\n", .{@tagName(ct)});
    }

    _vtable = cpu.vtableFromCoreType(ct);
    if (_vtable.init) |f| f();
    apple_accel.start(blas.config.accelerate);
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
pub const LogSink = struct {
    pub fn warn(_: LogSink, level: i32, comptime fmt: []const u8, args: anytype) void {
        blasWarn(level, fmt, args);
    }
};

/// A diagnostic, shown when `OPENBLAS_VERBOSE` is at least `level`. Only where
/// there is somewhere to print to.
fn blasWarn(level: i32, comptime fmt: []const u8, args: anytype) void {
    if (comptime sysinfo.target.os.tag == .freestanding) return;
    if (blas.config.verbose < level) return;
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
    blas.config = .{};
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
