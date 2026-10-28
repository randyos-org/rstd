//! Every level's portable kernels, which the vtable holds wherever a core has
//! none of its own; a core's kernels then go over them (`cpu.kernel.install`).
//! `kernel` picks a routine's kernel by element type: the vtable's where it
//! has a slot for the type, else the portable one built for it.

const std = @import("std");
const Vtable = @import("Vtable.zig");
const types = @import("types.zig");
const gemm_params = @import("blas3/gemm/params.zig");
const blas1 = @import("blas1/__root__.zig");
const blas2 = @import("blas2/__root__.zig");
const blas3 = @import("blas3/__root__.zig");
const bf16_kernel = @import("ext/bf16/kernel.zig");

const XDouble = types.XDouble;

const levels = .{ blas1, blas2, blas3 };

pub fn install(vt: *Vtable) void {
    // the generic target's blocking; `finish` sets the portable kernels' unrolls
    gemm_params.installCore(vt, .generic, .generic, .{});
    inline for (levels) |l| l.installGeneric(vt);
    finish(vt);
}

/// What follows from the kernels `vt` ends with; again after a core's:
/// the bfloat16-named slots with no kernel of their own point at the `f32`
/// kernels (`ext/bf16/kernel.zig`), and the gemm unrolls are those of the
/// portable kernels wherever those stayed (`blas3/gemm/kernel.zig`).
pub fn finish(vt: *Vtable) void {
    bf16_kernel.alias(vt);
    blas3.gemm.kernel.settle(vt);
}

// ---- type dispatch ---------------------------------------------------------

/// A kernel by the slot it fills, less the precision letter: `axpy_k` is
/// `saxpy_k` for `f32`; `iamax_k` is `isamax_k`, and complex `rot_k` is
/// `csrot_k`.
pub const Op = enum {
    copy_k,
    swap_k,
    scal_k,
    dot_k,
    dotu_k,
    dotc_k,
    nrm2_k,
    asum_k,
    sum_k,
    amax_k,
    amin_k,
    max_k,
    min_k,
    iamax_k,
    iamin_k,
    imax_k,
    imin_k,
    rot_k,
    rotm_k,
    axpy_k,
    axpyc_k,
    axpby_k,
    gemv_n,
    gemv_t,
    gemv_r,
    gemv_c,
    gemv_o,
    gemv_u,
    gemv_s,
    gemv_d,
    ger_k,
    geru_k,
    gerc_k,
    symv_U,
    symv_L,
    hemv_U,
    hemv_L,
    hemv_V,
    hemv_M,
    gemm_kernel,
    gemm_kernel_n,
    gemm_kernel_l,
    gemm_kernel_r,
    gemm_kernel_b,
    gemm_beta,
    gemm_incopy,
    gemm_itcopy,
    gemm_oncopy,
    gemm_otcopy,
    symm_iutcopy,
    symm_iltcopy,
    symm_outcopy,
    symm_oltcopy,
    hemm_iutcopy,
    hemm_iltcopy,
    hemm_outcopy,
    hemm_oltcopy,
    trmm_kernel_LN,
    trmm_kernel_LT,
    trmm_kernel_LR,
    trmm_kernel_LC,
    trmm_kernel_RN,
    trmm_kernel_RT,
    trmm_kernel_RR,
    trmm_kernel_RC,
    trmm_iunucopy,
    trmm_iunncopy,
    trmm_iutucopy,
    trmm_iutncopy,
    trmm_ilnucopy,
    trmm_ilnncopy,
    trmm_iltucopy,
    trmm_iltncopy,
    trmm_ounucopy,
    trmm_ounncopy,
    trmm_outucopy,
    trmm_outncopy,
    trmm_olnucopy,
    trmm_olnncopy,
    trmm_oltucopy,
    trmm_oltncopy,
    trsm_kernel_LN,
    trsm_kernel_LT,
    trsm_kernel_LR,
    trsm_kernel_LC,
    trsm_kernel_RN,
    trsm_kernel_RT,
    trsm_kernel_RR,
    trsm_kernel_RC,
    trsm_iunucopy,
    trsm_iunncopy,
    trsm_iutucopy,
    trsm_iutncopy,
    trsm_ilnucopy,
    trsm_ilnncopy,
    trsm_iltucopy,
    trsm_iltncopy,
    trsm_ounucopy,
    trsm_ounncopy,
    trsm_outucopy,
    trsm_outncopy,
    trsm_olnucopy,
    trsm_olnncopy,
    trsm_oltucopy,
    trsm_oltncopy,
    neg_tcopy,
};

/// The module whose `Portable(T)` builds `op`'s kernel.
fn Module(comptime op: Op) type {
    return switch (op) {
        .copy_k => blas1.copy,
        .swap_k => blas1.swap,
        .scal_k => blas1.scal,
        .dot_k, .dotu_k, .dotc_k => blas1.dot,
        .nrm2_k => blas1.nrm2,
        .asum_k => blas1.asum,
        .sum_k => blas1.sum,
        .amax_k, .amin_k, .max_k, .min_k, .iamax_k, .iamin_k, .imax_k, .imin_k => @import("blas1/extremum.zig"),
        .rot_k => blas1.rot,
        .rotm_k => blas1.rotm,
        .axpy_k, .axpyc_k => blas1.axpy,
        .axpby_k => blas1.axpby,
        .gemv_n, .gemv_t, .gemv_r, .gemv_c, .gemv_o, .gemv_u, .gemv_s, .gemv_d => blas2.gemv,
        .ger_k, .geru_k, .gerc_k => blas2.ger,
        .symv_U, .symv_L => blas2.symv,
        .hemv_U, .hemv_L, .hemv_V, .hemv_M => blas2.hemv,
        else => @compileError("level-3 kernels come from `blas3.gemm.kernel.find`"),
    };
}

/// The reference's precision letter for `T`, or null for a type with no
/// slots. `f64` comes before `XDouble`, which may be the same type.
pub fn letter(comptime T: type) ?[]const u8 {
    const F = types.Float(T);
    const c = types.isComplex(T);
    if (F == f32) return if (c) "c" else "s";
    if (F == f64) return if (c) "z" else "d";
    if (F == XDouble) return if (c) "x" else "q";
    return null;
}

/// Ops the reference has no extended-precision (`q`, `x`) slot for, which
/// run the portable kernel there: `axpby` exists for `s`, `d`, `c`, `z` only.
const no_extended_slot: std.EnumSet(Op) = .initMany(&.{.axpby_k});

/// The vtable slot `op` fills for `T`, or null for a type with no slots. A
/// type with slots that lacks this one is a compile error, unless
/// `no_extended_slot` says the reference has none.
pub fn slot(comptime T: type, comptime op: Op) ?[]const u8 {
    return comptime blk: {
        const l = letter(T) orelse break :blk null;
        const tag = @tagName(op);
        const name: []const u8 = switch (op) {
            .iamax_k, .iamin_k, .imax_k, .imin_k => "i" ++ l ++ tag[1..],
            .rot_k => if (types.isComplex(T)) l ++ letter(types.Float(T)).? ++ tag else l ++ tag,
            else => l ++ tag,
        };
        if (@hasField(Vtable, name)) break :blk name;
        const extended = std.mem.eql(u8, l, "q") or std.mem.eql(u8, l, "x");
        if (extended and no_extended_slot.contains(op)) break :blk null;
        @compileError("no vtable slot " ++ name ++ " for " ++ @typeName(T));
    };
}

/// The type of `op`'s kernel over `T`, any float type or its `Complex`.
pub fn Kernel(comptime T: type, comptime op: Op) type {
    if (@typeInfo(types.Float(T)) != .float)
        @compileError(@typeName(T) ++ " is not a float type or its Complex (bfloat16 and the mixed half routines have their own)");
    return PortableKernel(T, op);
}

/// Whether `op` is a level-3 kernel, which serves every routine over the
/// gemm tile (`blas3/gemm/kernel.zig`'s `sets`).
fn level3(comptime op: Op) bool {
    return @backingInt(op) >= @backingInt(Op.gemm_kernel);
}

fn PortableKernel(comptime T: type, comptime op: Op) type {
    if (level3(op)) return @TypeOf(blas3.gemm.kernel.find(T, @tagName(op)));
    return @TypeOf(@field(Module(op).Portable(T), @tagName(op)));
}

/// The portable kernel of `op` over `T`.
inline fn portableKernel(comptime T: type, comptime op: Op) PortableKernel(T, op) {
    if (comptime level3(op)) return blas3.gemm.kernel.find(T, @tagName(op));
    return @field(Module(op).Portable(T), @tagName(op));
}

/// `op`'s kernel for `T`: `vt`'s slot where it has one for `T`, else the
/// portable kernel over `T`.
pub inline fn kernel(comptime T: type, vt: *const Vtable, comptime op: Op) Kernel(T, op) {
    if (comptime slot(T, op)) |name| return @field(vt, name).?;
    return portableKernel(T, op);
}

/// `kernel` with the table from `getter` (`blas.vtable`, or `blas.table`
/// under a routine that fetched it), which a type without slots never calls.
pub inline fn fetch(comptime T: type, comptime op: Op, comptime getter: anytype) Kernel(T, op) {
    if (comptime slot(T, op)) |name| return @field(getter(), name).?;
    return portableKernel(T, op);
}

test "kernel takes the vtable's slot where there is one, else the portable kernel" {
    var vt: Vtable = .{};
    install(&vt);
    try std.testing.expect(kernel(f32, &vt, .axpy_k) == vt.saxpy_k.?);
    try std.testing.expect(kernel(std.math.Complex(f64), &vt, .rot_k) == vt.zdrot_k.?);
    try std.testing.expect(kernel(f64, &vt, .iamax_k) == vt.idamax_k.?);
    try std.testing.expectEqualStrings("csrot_k", slot(types.Complex32, .rot_k).?);
    // every type the tests use, by whether this build gives it slots
    inline for (.{ f16, f32, f64, f80, f128, XDouble }) |F| {
        inline for (.{ F, std.math.Complex(F) }) |T| {
            const got = kernel(T, &vt, .asum_k);
            const portable_k = @field(Module(.asum_k).Portable(T), "asum_k");
            if (comptime letter(T) == null) {
                try std.testing.expect(slot(T, .asum_k) == null);
                try std.testing.expect(got == portable_k);
                // Zig's own convention: only Zig calls it
                const E = types.Float(T);
                try std.testing.expect(@TypeOf(got) == *const fn (types.BlasLong, [*]E, types.BlasLong) E);
            } else {
                try std.testing.expect(got == @field(vt, slot(T, .asum_k).?).?);
            }
        }
    }
    // the documented gap: no extended axpby slot, so the portable kernel
    if (comptime XDouble != f64) try std.testing.expect(slot(XDouble, .axpby_k) == null);
}

test "every op has a slot for every type with slots, but the documented gaps" {
    @setEvalBranchQuota(20_000);
    const real_ops = [_]Op{ .copy_k, .swap_k, .scal_k, .dot_k, .nrm2_k, .asum_k, .sum_k, .amax_k, .amin_k, .max_k, .min_k, .iamax_k, .iamin_k, .imax_k, .imin_k, .rot_k, .rotm_k, .axpy_k, .axpby_k, .gemv_n, .gemv_t, .ger_k, .symv_U, .symv_L };
    const complex_ops = [_]Op{ .copy_k, .swap_k, .scal_k, .dotu_k, .dotc_k, .nrm2_k, .asum_k, .sum_k, .amax_k, .amin_k, .iamax_k, .iamin_k, .rot_k, .axpy_k, .axpyc_k, .axpby_k, .gemv_n, .gemv_t, .gemv_r, .gemv_c, .gemv_o, .gemv_u, .gemv_s, .gemv_d, .geru_k, .gerc_k, .symv_U, .symv_L, .hemv_U, .hemv_L, .hemv_V, .hemv_M };
    inline for (.{ f32, f64, XDouble }) |F| {
        inline for (.{ .{ F, real_ops }, .{ std.math.Complex(F), complex_ops } }) |p| {
            const T = p[0];
            inline for (p[1]) |op| {
                const gap = comptime no_extended_slot.contains(op) and F == XDouble and F != f64;
                try std.testing.expect((comptime slot(T, op) == null) == gap);
            }
        }
    }
}
