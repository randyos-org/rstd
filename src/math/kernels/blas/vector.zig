//! `@Vector` kernels for level-1 and level-2 slots, in the kernel ABI, `U`
//! native vectors a step; unit strides vectorize, others take a scalar loop,
//! and the special values are the portable kernels' (design.md). At `U = 4`
//! they mostly match or beat OpenBLAS's arm64 kernels on Apple's cores
//! (`tools/blasbench/results/results-m6.md`). `registry.Vector` places them.
//! The kernels are in `vector/` (`real.zig`, `complex.zig`, `half.zig`); here,
//! the slots they fill. Imports stay inside `linalg/blas/`, for the per-core
//! objects.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("types.zig");
const util = @import("util.zig");
const slots = @import("slots.zig");
const Vtable = @import("Vtable.zig");
const shared = @import("vector/util.zig");

const Slot = util.Slot;

pub const lanes = shared.lanes;
pub const splatTo = shared.splatTo;
pub const Real = @import("vector/real.zig").Real;
pub const Complex = @import("vector/complex.zig").Complex;
pub const Half = @import("vector/half.zig").Half;

/// Whether `slot`'s converted kernels give way to the portable or these on
/// every core, for one contract (`readme.md`, "One contract"): `nrm2`, the
/// extremes, `sdot` and `dsdot`, `.plain`'s single-precision `*_wide_k` and
/// the `peels` slots. Their assembly is
/// not converted (the indexes' `skipped`, kept in step).
pub fn routed(comptime slot: Slot) bool {
    return switch (slot) {
        // the extremes, `nrm2`, `sdot` and `dsdot`
        .samax_k,
        .damax_k,
        .camax_k,
        .zamax_k,
        .samin_k,
        .damin_k,
        .camin_k,
        .zamin_k,
        .smax_k,
        .dmax_k,
        .smin_k,
        .dmin_k,
        .isamax_k,
        .idamax_k,
        .icamax_k,
        .izamax_k,
        .isamin_k,
        .idamin_k,
        .icamin_k,
        .izamin_k,
        .ismax_k,
        .idmax_k,
        .ismin_k,
        .idmin_k,
        .snrm2_k,
        .dnrm2_k,
        .cnrm2_k,
        .znrm2_k,
        .sdot_k,
        .dsdot_k,
        // `.plain`'s single-precision sums in `f64`, which no converted kernel fills (issue #584)
        .sasum_wide_k,
        .casum_wide_k,
        .ssum_wide_k,
        .csum_wide_k,
        .cdotu_wide_k,
        .cdotc_wide_k,
        => true,
        else => peels(slot),
    };
}

/// The x86_64 reductions whose converted kernels peel elements until an
/// operand is aligned, so their bits follow its address: `asum_sse.S`,
/// `asum_sse2.S`, `zasum_sse*.S`, `zsum_sse*.S`, `dot_sse2.S`,
/// `ddot_bulldozer.S`, `zdot_sse*.S`, `dgemv_n.S` and their `atom` and
/// `bulldozer` kin. Their C kernels on other cores stay.
pub fn peels(comptime slot: Slot) bool {
    if (sysinfo.target.cpu.arch != .x86_64) return false;
    return std.mem.findScalar(Slot, &peeling, slot) != null;
}

/// The slots `peels` names.
pub const peeling = [_]Slot{ .sasum_k, .dasum_k, .casum_k, .zasum_k, .csum_k, .zsum_k, .ddot_k, .cdotu_k, .cdotc_k, .zdotu_k, .zdotc_k, .dgemv_n };

/// A function of `Real` or `Complex`, by its name there.
const Fn = enum {
    axpy,
    axpyc,
    dot,
    dsdot,
    dotu,
    dotc,
    scal,
    nrm2,
    asum,
    sum,
    asum_wide,
    sum_wide,
    dotu_wide,
    dotc_wide,
    copy,
    swap,
    rot,
    gemv_n,
    gemv_t,
    gemv_r,
    gemv_c,
    gemv_o,
    gemv_u,
    gemv_s,
    gemv_d,
    amax,
    iamax,
    amin,
    iamin,
    max,
    imax,
    min,
    imin,

    /// Whether it is one of `.plain`'s single-precision sums in `f64` (issue #584).
    fn wide(f: Fn) bool {
        return switch (f) {
            .asum_wide, .sum_wide, .dotu_wide, .dotc_wide => true,
            else => false,
        };
    }

    /// Its op; null for `dsdot`, which has one slot (`dsdot_k`), not an op.
    fn op(f: Fn) ?slots.Op {
        return switch (f) {
            .axpy => .axpy_k,
            .axpyc => .axpyc_k,
            .dot => .dot_k,
            .dsdot => null,
            .dotu => .dotu_k,
            .dotc => .dotc_k,
            .scal => .scal_k,
            .nrm2 => .nrm2_k,
            .asum => .asum_k,
            .sum => .sum_k,
            .asum_wide => .asum_wide_k,
            .sum_wide => .sum_wide_k,
            .dotu_wide => .dotu_wide_k,
            .dotc_wide => .dotc_wide_k,
            .copy => .copy_k,
            .swap => .swap_k,
            .rot => .rot_k,
            .gemv_n => .gemv_n,
            .gemv_t => .gemv_t,
            .gemv_r => .gemv_r,
            .gemv_c => .gemv_c,
            .gemv_o => .gemv_o,
            .gemv_u => .gemv_u,
            .gemv_s => .gemv_s,
            .gemv_d => .gemv_d,
            .amax => .amax_k,
            .iamax => .iamax_k,
            .amin => .amin_k,
            .iamin => .iamin_k,
            .max => .max_k,
            .imax => .imax_k,
            .min => .min_k,
            .imin => .imin_k,
        };
    }
};

/// The functions of `Real`, each filling the slot of its op (`slotOf`).
const real_fns = [_]Fn{ .axpy, .dot, .dsdot, .scal, .nrm2, .asum, .sum, .asum_wide, .sum_wide, .copy, .swap, .rot, .gemv_n, .gemv_t, .amax, .iamax, .amin, .iamin, .max, .imax, .min, .imin };

/// The functions of `Complex`, likewise.
const complex_fns = [_]Fn{ .axpy, .axpyc, .dotu, .dotc, .scal, .nrm2, .asum, .sum, .asum_wide, .sum_wide, .dotu_wide, .dotc_wide, .copy, .swap, .amax, .iamax, .amin, .iamin, .rot, .gemv_n, .gemv_t, .gemv_r, .gemv_c, .gemv_o, .gemv_u, .gemv_s, .gemv_d };

/// The slot function `f` over `T` fills: its op at `T`'s letter, or null for
/// `dsdot` past `f32` and for the `*_wide` functions past single precision.
fn slotOf(comptime T: type, comptime f: Fn) ?Slot {
    if (f.wide() and types.Float(T) != f32) return null;
    const op = f.op() orelse return if (T == f32) .dsdot_k else null;
    return slots.slotTag(.from(slots.letter(T).?), op);
}

/// A kernel: the slot it fills, and the function of `Real` or `Complex` over
/// `F` that is it.
const Kernel = struct { slot: Slot, F: type, complex: bool, f: Fn };

/// Every function here as a kernel, whatever its type's lanes do with
/// denormals.
pub const every: []const Kernel = blk: {
    // every function at each type, past the default quota
    @setEvalBranchQuota(10_000);
    var r: []const Kernel = &.{};
    for ([_]type{ f32, f64 }) |F| for (real_fns) |f| {
        if (slotOf(F, f)) |s| r = r ++ &[_]Kernel{.{ .slot = s, .F = F, .complex = false, .f = f }};
    };
    for ([_]type{ f32, f64 }) |F| for (complex_fns) |f| {
        if (slotOf(std.math.Complex(F), f)) |s| r = r ++ &[_]Kernel{.{ .slot = s, .F = F, .complex = true, .f = f }};
    };
    break :blk r;
};

/// Every kernel here: none in a type whose lanes flush denormals
/// (`util.vectorIeee`), whose slots keep the portable kernels.
pub const kernels: []const Kernel = blk: {
    var r: []const Kernel = &.{};
    for (every) |k| {
        if (util.vectorIeee(k.F)) r = r ++ &[_]Kernel{k};
    }
    break :blk r;
};

/// The kernel of `k`, `U` vectors a step.
pub fn kernel(comptime k: Kernel, comptime U: comptime_int) util.SlotFn(k.slot) {
    return &@field(if (k.complex) Complex(k.F, U) else Real(k.F, U), @tagName(k.f));
}

/// The symbol of `slot`'s kernel here for `core`, apart from a converted one's.
pub fn symbol(comptime slot: Slot, comptime core: anytype) []const u8 {
    return util.kernelSymbol(slot, core) ++ "_vec";
}

/// Points every slot these kernels fill at this build's own, four vectors a
/// step, `ported`'s too, for the tests.
pub fn installAll(vt: *Vtable) void {
    inline for (kernels) |k| @field(vt, @tagName(k.slot)) = kernel(k, 4);
    inline for (ported.slots) |s| {
        if (comptime ported.builds(s)) @field(vt, @tagName(s)) = ported.entry(s);
    }
}

/// The kernels here no `Registry` routes, outside `real_slots` and
/// `complex_slots`: an architecture's hand ports name them (`entry`, as a
/// port's module) where upstream has intrinsics for them.
pub const ported = struct {
    pub const slots = [_]Slot{
        .saxpby_k, .daxpby_k, .caxpby_k, .zaxpby_k,
        .srotm_k,  .drotm_k,  .ssymv_U,  .ssymv_L,
        .dsymv_U,  .dsymv_L,  .csymv_U,  .csymv_L,
        .zsymv_U,  .zsymv_L,  .chemv_U,  .chemv_L,
        .chemv_V,  .chemv_M,  .zhemv_U,  .zhemv_L,
        .zhemv_V,  .zhemv_M,  .sbgemv_n, .sbgemv_t,
        .shgemv_n, .shgemv_t, .bgemv_n,  .bgemv_t,
    };

    /// `slot`'s kernel, four vectors a step.
    pub fn entry(comptime slot: Slot) util.SlotFn(slot) {
        const R = Real(f32, 4);
        const D = Real(f64, 4);
        const C = Complex(f32, 4);
        const Z = Complex(f64, 4);
        return switch (slot) {
            .saxpby_k => &R.axpby,
            .daxpby_k => &D.axpby,
            .caxpby_k => &C.axpby,
            .zaxpby_k => &Z.axpby,
            .srotm_k => &R.rotm,
            .drotm_k => &D.rotm,
            .ssymv_U => &R.symv_U,
            .ssymv_L => &R.symv_L,
            .dsymv_U => &D.symv_U,
            .dsymv_L => &D.symv_L,
            .csymv_U => &C.symv_U,
            .csymv_L => &C.symv_L,
            .zsymv_U => &Z.symv_U,
            .zsymv_L => &Z.symv_L,
            .chemv_U => &C.hemv_U,
            .chemv_L => &C.hemv_L,
            .chemv_V => &C.hemv_V,
            .chemv_M => &C.hemv_M,
            .zhemv_U => &Z.hemv_U,
            .zhemv_L => &Z.hemv_L,
            .zhemv_V => &Z.hemv_V,
            .zhemv_M => &Z.hemv_M,
            .sbgemv_n => &Half(types.BFloat16, f32).gemv_n,
            .sbgemv_t => &Half(types.BFloat16, f32).gemv_t,
            .shgemv_n => &Half(types.HF16, f32).gemv_n,
            .shgemv_t => &Half(types.HF16, f32).gemv_t,
            .bgemv_n => &Half(types.BFloat16, types.BFloat16).gemv_n,
            .bgemv_t => &Half(types.BFloat16, types.BFloat16).gemv_t,
            else => @compileError(@tagName(slot) ++ " has no kernel in vector.ported"),
        };
    }

    /// Whether this build makes `slot`'s kernel for `installAll`: none over
    /// `f32` whose lanes flush (`util.vectorIeee`), and the half-precision
    /// ones only where a port uses them or the native tests run, since
    /// `f16` vectors stop LLVM on some targets.
    pub fn builds(comptime slot: Slot) bool {
        const name = @tagName(slot);
        const half = for ([_][]const u8{ "sbgemv", "shgemv", "bgemv" }) |f| {
            if (std.mem.startsWith(u8, name, f)) break true;
        } else false;
        if (half) return switch (sysinfo.target.cpu.arch) {
            .x86_64, .aarch64, .riscv64 => true,
            else => false,
        };
        const single = name[0] == 's' or name[0] == 'c';
        return !single or util.vectorIeee(f32);
    }
};
