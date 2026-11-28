//! Apple Accelerate's CBLAS, bound at run time through `dll.zig` with no SDK,
//! and the routing to it (design.md, Apple Accelerate). Sources:
//! <https://developer.apple.com/documentation/macos-release-notes/macos-13_3-release-notes>,
//! <https://developer.apple.com/documentation/accelerate/blas-library>,
//! <https://developer.apple.com/documentation/accelerate/blassetthreading(_:)>.
//!
//! Here, the routing and the call; in `apple_accel/`, the table and its
//! loading (`binding.zig`), a call's parts (`split.zig`), `quirks.zig`, the
//! speed tables (`speed.zig`), and
//! the routines by level (`level1.zig` to `level3.zig`, `lapack.zig`).

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("types.zig");
const blas = @import("__root__.zig");
const binding = @import("apple_accel/binding.zig");
const routines = @import("apple_accel/routines.zig");
const parts = @import("apple_accel/split.zig");
const quirks = @import("apple_accel/quirks.zig");
const speed = @import("apple_accel/speed.zig");
const level1 = @import("apple_accel/level1.zig");
const level2 = @import("apple_accel/level2.zig");
const level3 = @import("apple_accel/level3.zig");
const lapack = @import("apple_accel/lapack.zig");

const Complex32 = types.Complex32;
const Complex64 = types.Complex64;

pub const available = binding.available;
pub const x86_64 = binding.x86_64;
pub const potrf_scratches_lower = binding.potrf_scratches_lower;
pub const Order = binding.Order;
pub const Transpose = binding.Transpose;
pub const Uplo = binding.Uplo;
pub const Diag = binding.Diag;
pub const Side = binding.Side;
pub const Threading = binding.Threading;
pub const Table = binding.Table;
pub const ilp64Names = binding.ilp64Names;
pub const lp64Names = binding.lp64Names;
pub const fortran = binding.fortran;
pub const ThreadApi = binding.ThreadApi;
pub const Blas = binding.Blas;
pub const AppleAccelerate = binding.AppleAccelerate;
pub const Options = binding.Options;
pub const bindBlas = binding.bindBlas;
pub const load = binding.load;
pub const EnvError = binding.EnvError;
pub const setMaxThreadsEnv = binding.setMaxThreadsEnv;
pub const Symbol = routines.Symbol;

pub const vector_quantum = parts.vector_quantum;
pub const matrix_quantum = parts.matrix_quantum;
pub const L1 = speed.L1;
pub const sent = speed.sent;
pub const L23 = speed.L23;
pub const routes = speed.routes;

pub const copy = level1.copy;
pub const swap = level1.swap;
pub const scal = level1.scal;
pub const dot = level1.dot;
pub const dsdot = level1.dsdot;
pub const nrm2 = level1.nrm2;
pub const asum = level1.asum;
pub const rot = level1.rot;
pub const rotm = level1.rotm;
pub const axpy = level1.axpy;
pub const axpby = level1.axpby;

pub const gemv = level2.gemv;
pub const gbmv = level2.gbmv;
pub const ger = level2.ger;
pub const symv = level2.symv;
pub const sbmv = level2.sbmv;
pub const TriKind = level2.TriKind;
pub const tri = level2.tri;
pub const syr = level2.syr;
pub const syr2 = level2.syr2;

pub const gemm = level3.gemm;
pub const symm = level3.symm;
pub const syrk = level3.syrk;
pub const trmm = level3.trmm;

pub const potrf = lapack.potrf;
pub const getrf = lapack.getrf;
pub const getrs = lapack.getrs;

// ---- routing -----------------------------------------------------------------
//
// Called after the public routines' own checks, so the errors are rstd's;
// false (or null) where rstd's kernels run (design.md, Apple Accelerate).

/// Whether the public routines over `T` can go to Accelerate.
pub fn covers(comptime T: type) bool {
    return available and (T == f32 or T == f64 or T == Complex32 or T == Complex64);
}

/// Whether this x86_64 process runs under Rosetta 2 on an Apple Silicon Mac,
/// read once by `start` (`sysctl.proc_translated`): its quirks (design.md)
/// keep `sdot`, `dsdot`, `snrm2`, `scnrm2` and `gemm` with a non-finite `A`
/// on rstd's kernels, which native x86_64 sends to Accelerate (issue #388).
var translated = false;

/// Whether Rosetta's quirks apply: an x86_64 build under Rosetta 2.
pub fn rosetta() bool {
    return x86_64 and translated;
}

fn readTranslated() bool {
    if (comptime !x86_64) return false;
    var v: c_int = 0;
    var len: usize = @sizeOf(c_int);
    if (std.c.sysctlbyname("sysctl.proc_translated", &v, &len, null, 0) != 0) return false;
    return v == 1;
}

/// The framework, loaded once by `start` and never unloaded.
pub var loaded: ?AppleAccelerate = null;
var tried = false;
/// Whether the public routines call it.
var on = std.atomic.Value(bool).init(false);
/// Whether this thread's calls are single-threaded (`singleThreadHere`).
threadlocal var thread_set = false;
/// Calls made into Accelerate, counted in tests only.
pub var calls = std.atomic.Value(usize).init(0);
/// Calls on a thread not single-threaded (macOS 15+), counted in tests only.
pub var unpinned = std.atomic.Value(usize).init(0);

/// `blas.init`'s part: loads Accelerate the first time `enable` holds (it
/// `setenv`s, so only here), then routes to it if it loaded.
pub fn start(enable: bool) void {
    if (comptime !available) return;
    if (enable and !tried) {
        tried = true;
        translated = readTranslated();
        loaded = load(.{});
    }
    on.store(enable and loaded != null, .release);
}

/// Routes to Accelerate or stops; true if it is in use after. Never loads
/// (only `start` does) nor unloads.
pub fn route(enable: bool) bool {
    const now = enable and loaded != null;
    on.store(now, .release);
    return now;
}

/// Whether the public routines call Accelerate now.
pub fn routing() bool {
    return on.load(.acquire);
}

/// The table for `T`'s routines, or null where rstd runs them.
pub inline fn use(comptime T: type) ?*const Blas {
    if (comptime !covers(T)) return null;
    _ = blas.vtable();
    if (!on.load(.acquire)) return null;
    return &loaded.?.blas;
}

/// Makes the calling thread's calls single-threaded, on its first.
fn pin() void {
    if (thread_set) return;
    thread_set = true;
    if (loaded) |*acc| acc.singleThreadHere();
}

/// What `call` returns for `sym`: its C result, an index as `isize`.
fn Ret(comptime sym: Symbol) type {
    const R = @typeInfo(@typeInfo(@typeInfo(@FieldType(Table(i64), @tagName(sym))).optional.child).pointer.child).@"fn".return_type.?;
    return if (R == i64) isize else R;
}

/// Calls `sym` with `args`, every `isize` narrowed to the table's integer;
/// null if the table lacks it or, on LP64, any argument or `extents` entry
/// (an operand's span in elements) is past `c_int`. Pins the calling thread
/// (`pin`).
pub fn call(comptime sym: Symbol, b: *const Blas, extents: []const isize, args: anytype) ?Ret(sym) {
    switch (b.*) {
        inline else => |*bound, tag| {
            const f = @field(bound.fns, @tagName(sym)) orelse return null;
            const I = if (tag == .ilp64) i64 else c_int;
            if (I != i64) for (extents) |e| if (e > std.math.maxInt(I)) return null;
            var conv: std.meta.ArgsTuple(@typeInfo(@TypeOf(f)).pointer.child) = undefined;
            inline for (0..args.len) |i| {
                conv[i] = if (@TypeOf(args[i]) == isize) std.math.cast(I, args[i]) orelse return null else args[i];
            }
            pin();
            if (sysinfo.is_test) {
                _ = calls.fetchAdd(1, .monotonic);
                if (loaded) |*acc| if (acc.threading()) |t| if (t != .single_threaded) {
                    _ = unpinned.fetchAdd(1, .monotonic);
                };
            }
            const r = @call(.auto, f, conv);
            return if (@TypeOf(r) == I) @intCast(r) else r;
        },
    }
}

/// The span of `n` elements at stride `inc`.
pub fn span(n: isize, inc: isize) isize {
    return (n - 1) * @as(isize, @intCast(@abs(inc))) + 1;
}

/// The span of a `rows` by `cols` matrix at leading dimension `ld`.
pub fn mat(rows: isize, cols: isize, ld: isize) isize {
    return ld * (cols - 1) + rows;
}
