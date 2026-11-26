//! x86_64 `ger` kernels: the `@Vector` port of `kernel/x86_64/sger.c` and
//! `dger.c` with Sandy Bridge's blocks (`sger_microk_sandy-2.c`,
//! `dger_microk_sandy-2.c`): per column, `a += (alpha*y[j])*x` unfused over
//! rows in sixteens, the rest by the core's own `axpy` kernel, as `AXPYU_K`.
//! Each element is one product and one sum, so the vector width leaves the
//! bits as they are. Tests live in `ger/__root__.zig`.

const coretypes = @import("../../../cpu/coretypes.zig");
const registry = @import("../../../registry.zig");
const types = @import("../../../types.zig");
const util = @import("../../../util.zig");

const BlasLong = types.BlasLong;
const Slot = util.Slot;

pub const Core = coretypes.X86CoreType;

/// `?ger.c` of `core`, with its `?ger_kernel_16`.
fn Ger(comptime F: type, comptime core: Core, comptime axpy: Slot) type {
    // the block's 128-bit registers, four at a time
    const V = @Vector(64 / @sizeOf(F), F);
    const l = @typeInfo(V).vector.len;
    return struct {
        fn k16(n: usize, x: [*]const F, a: [*]F, y0: F) void {
            const yv: V = @splat(y0);
            var i: usize = 0;
            while (i < n) : (i += l) util.store(V, a + i, util.load(V, a + i) + util.load(V, x + i) * yv);
        }

        pub fn k(m: BlasLong, n: BlasLong, _: BlasLong, alpha: F, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, a_in: [*]F, lda: BlasLong, buffer: [*]F) callconv(.c) i32 {
            if (m <= 0 or n <= 0) return 0;
            const mu: usize = @intCast(m);
            var xx = x;
            if (incx != 1) {
                xx = buffer;
                for (0..mu) |i| xx[i] = util.at(F, x, @as(BlasLong, @intCast(i)) * incx).*;
            }
            const m1 = mu & ~@as(usize, 15);
            const step: usize = @intCast(lda);
            const tail = registry.externKernel(axpy, core);
            var a = a_in;
            var yp = y;
            for (0..@intCast(n)) |_| {
                const y0 = alpha * yp[0];
                if (m1 > 0) k16(m1, xx, a, y0);
                if (mu > m1) _ = tail(@intCast(mu - m1), 0, 0, y0, xx + m1, 1, a + m1, 1, xx, 0);
                a += step;
                yp = @ptrCast(util.at(F, yp, incy));
            }
            return 0;
        }
    };
}

/// Sandy Bridge's kernels.
const Sandy = struct {
    pub fn entry(comptime slot: Slot) util.SlotFn(slot) {
        return switch (slot) {
            .sger_k => &Ger(f32, .sandybridge, .saxpy_k).k,
            .dger_k => &Ger(f64, .sandybridge, .daxpy_k).k,
            else => @compileError(@tagName(slot) ++ " is not a real ger slot"),
        };
    }
};

/// The ports and where they serve (`KERNEL.<CORE>`'s `?GERKERNEL`).
pub const ports = [_]registry.Port(Core){
    .{ .module = Sandy, .cores = &.{.sandybridge}, .slots = &.{ .sger_k, .dger_k }, .written = .vector },
};

const P = registry.Ports(Core, "x86_64", &ports);

pub const hasSlot = P.hasSlot;
pub const has = P.has;
pub const exportCore = P.exportCore;
pub const install = P.install;
