//! Apple Accelerate's vDSP DFTs (`vDSP_DFT_zop`, `vDSP_DFT_zrop`, macOS 10.7+):
//! `Plan` and `RealPlan` hand vDSP's lengths, `f 2^k` with `f` in 1, 3, 5
//! and 15, to it on macOS where it ran faster than rstd's own (`faster`), as
//! the BLAS hands its calls to Accelerate (issues #155, #856).
//! `useAccelerate(false)` before planning keeps rstd's own.
//! vDSP takes split complex data; one call per transform, single-threaded,
//! so a batch's bits depend on neither the worker count nor the order.

const std = @import("std");
const sysinfo = @import("builtin");
const dll = @import("../../../dll.zig");
const cmplx = @import("cmplx.zig");
const pffft = @import("pffft/__root__.zig");

const C = cmplx.C;

/// Where vDSP can be loaded.
pub const available = sysinfo.target.os.tag == .macos;

/// `vDSP_DFT_Direction`.
const Dir = enum(c_int) { forward = 1, inverse = -1 };

const Create = ?*const fn (?*anyopaque, usize, Dir) callconv(.c) ?*anyopaque;
fn Execute(comptime F: type) type {
    return ?*const fn (*anyopaque, [*]const F, [*]const F, [*]F, [*]F) callconv(.c) void;
}
const Destroy = ?*const fn (?*anyopaque) callconv(.c) void;

/// The vDSP routines used, by their symbols (`vDSP.h`).
const Table = struct {
    vDSP_DFT_zop_CreateSetup: Create,
    vDSP_DFT_zop_CreateSetupD: Create,
    vDSP_DFT_zrop_CreateSetup: Create,
    vDSP_DFT_zrop_CreateSetupD: Create,
    vDSP_DFT_Execute: Execute(f32),
    vDSP_DFT_ExecuteD: Execute(f64),
    vDSP_DFT_DestroySetup: Destroy,
    vDSP_DFT_DestroySetupD: Destroy,
};

const State = enum(u8) { unloaded, loading, loaded, failed };
var state = std.atomic.Value(State).init(.unloaded);
var table: Table = undefined;
/// Whether new plans use vDSP.
var on = std.atomic.Value(bool).init(true);
/// Whether new plans use it at every length it takes, not only `faster` ones.
var every = std.atomic.Value(bool).init(false);

/// Turns the use of vDSP by plans made after this on or off (macOS, `f32`
/// and `f64`); true if new plans use it. Plans keep what they were made with.
pub fn useAccelerate(enable: bool) bool {
    on.store(enable, .release);
    return enable and loaded() != null;
}

/// Makes plans made after this take vDSP at every length it takes, `faster`
/// or not (`enable`), or only where it is (the default): for the tests and
/// the bench's `fft` section, which compare them all.
pub fn useAccelerateEverywhere(enable: bool) void {
    every.store(enable, .release);
}

/// Whether vDSP ran length `n` (one of `validLength`'s) of `Plan(T)`, or of
/// `RealPlan(T)` if `real`, at least 0.97 times rstd's own in two runs of
/// the bench's `fft` section. x86_64 (the Intel Mac, #155): everywhere but
/// complex `f32` at pffft's sizes. aarch64 (the M6, #856): complex `f32`
/// from 3840 and where pffft lacks `n`, complex `f64` at powers of 2, real
/// from 32, `f64` only at powers of 2 up to 16384. Not complex `f32` past
/// 65536 but at powers of 2 on aarch64, where vDSP's error from an `f64`
/// transform was 18 to 24 eps of the largest entry, pffft's 2 to 4 (#870).
pub fn faster(comptime T: type, comptime real: bool, n: usize) bool {
    const pow2 = std.math.isPowerOfTwo(n);
    if (comptime sysinfo.target.cpu.arch == .aarch64) return switch (T) {
        f32 => if (real) n >= 32 else (n >= 3840 and (pow2 or n <= 65536)) or !pffft.validSize(n, .complex),
        f64 => if (real) pow2 and n >= 32 and n <= 16384 else pow2,
        else => false,
    };
    return real or T != f32 or !pffft.validSize(n, .complex);
}

/// The bound table, loading Accelerate the first time; null off macOS or
/// where a symbol is missing.
fn loaded() ?*const Table {
    if (comptime !available) return null;
    while (true) switch (state.load(.acquire)) {
        .loaded => return &table,
        .failed => return null,
        .loading => std.atomic.spinLoopHint(),
        .unloaded => if (state.cmpxchgStrong(.unloaded, .loading, .acquire, .monotonic) == null) {
            // kept open for the process: the setups point into it
            var fw = dll.openDarwinFramework("Accelerate") catch {
                state.store(.failed, .release);
                return null;
            };
            const b = dll.bind(Table, &fw, dll.exactName);
            table = b.fns;
            state.store(if (b.count() == @typeInfo(Table).@"struct".field_names.len) .loaded else .failed, .release);
        },
    };
}

/// One length's forward and inverse setups for `F`, complex or real
/// (`zrop`), or null where vDSP is off, absent or rejects the length.
pub fn Setup(comptime F: type, comptime real: bool) type {
    if (F != f32 and F != f64) @compileError("vDSP: f32 and f64 only");
    return struct {
        const Self = @This();
        fwd: *anyopaque,
        inv: *anyopaque,
        exec: *const fn (*anyopaque, [*]const F, [*]const F, [*]F, [*]F) callconv(.c) void,
        destroy: *const fn (?*anyopaque) callconv(.c) void,

        pub fn init(n: usize) ?Self {
            if (!on.load(.acquire) or !validLength(n, real)) return null;
            const t = loaded() orelse return null;
            const create = (if (F == f32)
                (if (real) t.vDSP_DFT_zrop_CreateSetup else t.vDSP_DFT_zop_CreateSetup)
            else if (real) t.vDSP_DFT_zrop_CreateSetupD else t.vDSP_DFT_zop_CreateSetupD).?;
            const destroy = (if (F == f32) t.vDSP_DFT_DestroySetup else t.vDSP_DFT_DestroySetupD).?;
            const f = create(null, n, .forward) orelse return null;
            // the inverse shares the forward's tables
            const i = create(f, n, .inverse) orelse {
                destroy(f);
                return null;
            };
            return .{ .fwd = f, .inv = i, .exec = (if (F == f32) t.vDSP_DFT_Execute else t.vDSP_DFT_ExecuteD).?, .destroy = destroy };
        }

        pub fn deinit(s: *Self) void {
            s.destroy(s.inv);
            s.destroy(s.fwd);
            s.* = undefined;
        }

        /// The split input `ir`, `ii` to `or_`, `oi`, unscaled.
        pub fn run(s: *const Self, fwd: bool, ir: []const F, ii: []const F, or_: []F, oi: []F) void {
            s.exec(if (fwd) s.fwd else s.inv, ir.ptr, ii.ptr, or_.ptr, oi.ptr);
        }
    };
}

/// The lengths planned on vDSP: `f 2^k`, `f` in 1, 3, 5 and 15 (`vDSP.h`), `k` at
/// least 3 complex and 4 real; vDSP's setup may still refuse one (`f64` takes
/// fewer `15 2^k` than `f32`), which `Setup.init` reports as null.
pub fn validLength(n: usize, real: bool) bool {
    var m = n;
    var k: usize = 0;
    while (m != 0 and m % 2 == 0) : (k += 1) m /= 2;
    return (m == 1 or m == 3 or m == 5 or m == 15) and k >= @as(usize, if (real) 4 else 3);
}

/// `Plan`'s and `RealPlan`'s vDSP setup for `T`: `Setup` for `f32` and
/// `f64` at the lengths it is `faster` at, a type no plan holds for the others.
pub fn Vdsp(comptime T: type, comptime real: bool) type {
    if (T != f32 and T != f64) return struct {
        const Self = @This();
        pub fn init(_: usize) ?Self {
            return null;
        }
        pub fn deinit(_: *Self, _: std.mem.Allocator) void {}
        pub fn run(_: *const Self, _: []C(T), _: []C(T), _: bool, _: T) void {}
        pub fn runReal(_: *const Self, _: []T, _: []T, _: bool, _: T) void {}
    };
    return struct {
        const Self = @This();
        setup: Setup(T, real),

        pub fn init(n: usize) ?Self {
            if (!every.load(.acquire) and !faster(T, real, n)) return null;
            return .{ .setup = Setup(T, real).init(n) orelse return null };
        }

        pub fn deinit(self: *Self, _: std.mem.Allocator) void {
            self.setup.deinit();
        }

        /// `c` transformed in place and scaled by `f`; `rest` holds `2 n`.
        pub fn run(self: *const Self, c: []C(T), rest: []C(T), fwd: bool, f: T) void {
            const n = c.len;
            const fl: [*]T = @ptrCast(rest.ptr);
            const ir = fl[0..n];
            const ii = fl[n .. 2 * n];
            const or_ = fl[2 * n .. 3 * n];
            const oi = fl[3 * n .. 4 * n];
            for (c, ir, ii) |x, *r, *i| {
                r.* = x.re;
                i.* = x.im;
            }
            self.setup.run(fwd, ir, ii, or_, oi);
            for (c, or_, oi) |*x, r, i| x.* = .{ .re = r * f, .im = i * f };
        }

        /// `buf`, `n` reals in FFTPACK's halfcomplex order (`RealPlan`'s), transformed
        /// and scaled by `f`; `rest` holds `2 n`. vDSP takes the reals as even and
        /// odd halves and packs the half spectrum with the Nyquist value in element
        /// 0's imaginary part, twice the transform forward and `n` times back.
        pub fn runReal(self: *const Self, buf: []T, rest: []T, fwd: bool, f: T) void {
            const n = buf.len;
            const h = n / 2;
            const ir = rest[0..h];
            const ii = rest[h..n];
            const or_ = rest[n .. n + h];
            const oi = rest[n + h .. 2 * n];
            if (fwd) {
                for (ir, ii, 0..) |*r, *i, j| {
                    r.* = buf[2 * j];
                    i.* = buf[2 * j + 1];
                }
                self.setup.run(true, ir, ii, or_, oi);
                const g = f / 2;
                buf[0] = or_[0] * g;
                buf[n - 1] = oi[0] * g;
                for (1..h) |k| {
                    buf[2 * k - 1] = or_[k] * g;
                    buf[2 * k] = oi[k] * g;
                }
            } else {
                ir[0] = buf[0];
                ii[0] = buf[n - 1];
                for (1..h) |k| {
                    ir[k] = buf[2 * k - 1];
                    ii[k] = buf[2 * k];
                }
                self.setup.run(false, ir, ii, or_, oi);
                for (or_, oi, 0..) |r, i, j| {
                    buf[2 * j] = r * f;
                    buf[2 * j + 1] = i * f;
                }
            }
        }
    };
}
