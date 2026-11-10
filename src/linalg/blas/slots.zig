//! Vtable slot names by enum: the gemms, their blocking parameters and
//! kernel slots, and the half-precision gemvs. A leaf, as the kernel objects
//! compile it (`ext/bf16/kernel.zig`); `portable.zig` re-exports it.

const types = @import("types.zig");
const util = @import("util.zig");
const Vtable = @import("Vtable.zig");

const XDouble = types.XDouble;

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

/// A gemm with blocking parameters, by `blas3/gemm/params.zig`'s `Table`
/// field: a precision, a half-precision gemm or a 3M one.
pub const Gemm = enum {
    s,
    d,
    q,
    c,
    z,
    x,
    sh,
    b,
    sb,
    c3m,
    z3m,
    x3m,

    /// The 3M gemm of complex `g`.
    pub fn threeM(comptime g: Gemm) Gemm {
        return switch (g) {
            .c => .c3m,
            .z => .z3m,
            .x => .x3m,
            else => @compileError(@tagName(g) ++ " has no 3M gemm"),
        };
    }

    /// `T`'s gemm, its 3M one for `three`.
    pub fn of(comptime T: type, comptime three: bool) Gemm {
        return comptime @field(Gemm, letter(T).? ++ if (three) "3m" else "");
    }

    /// Its vtable fields' prefix: `sgemm`, `cgemm3m`, `sbgemm`.
    fn prefix(comptime g: Gemm) []const u8 {
        return comptime switch (g) {
            .c3m, .z3m, .x3m => @tagName(g)[0..1] ++ "gemm3m",
            else => @tagName(g) ++ "gemm",
        };
    }
};

/// A gemm's blocking parameter: `.p` is `sgemm_p` for `.s`, `cgemm3m_p`
/// for `.c3m`.
pub const Param = enum { p, q, r, unroll_m, unroll_n, unroll_mn, align_k };

/// `paramPtr`'s `Slot`, for a kernel object's symbol.
pub fn paramTag(comptime g: Gemm, comptime param: Param) util.Slot {
    return comptime @field(util.Slot, g.prefix() ++ "_" ++ @tagName(param));
}

/// The slot of precision `g`'s level-3 kernel `name`, as a kernel set
/// names it (the slot less its letter: `gemm_kernel` is `dgemm_kernel` for
/// `.d`), or null if the vtable has none.
pub fn kernelSlot(comptime g: Gemm, comptime name: []const u8) ?util.Slot {
    const field = @tagName(g) ++ name;
    return comptime if (@hasField(Vtable, field)) @field(util.Slot, field) else null;
}

/// `vt`'s field for `s`, by pointer; `vt` points to a `Vtable`.
pub inline fn fieldPtr(vt: anytype, comptime s: util.Slot) @TypeOf(&@field(vt, @tagName(s))) {
    return &@field(vt, @tagName(s));
}

/// `vt`'s parameter `param` of `g`, by pointer; `vt` points to a `Vtable`.
pub inline fn paramPtr(vt: anytype, comptime g: Gemm, comptime param: Param) @TypeOf(&@field(vt, g.prefix() ++ "_" ++ @tagName(param))) {
    return &@field(vt, g.prefix() ++ "_" ++ @tagName(param));
}

/// A gemm's kernel slot, less its gemm's prefix: `.kernel` is
/// `sbgemm_kernel` for `.sb`, `.small_kernel_b0_nt` `dgemm_small_kernel_b0_nt`
/// for `.d`.
pub const GemmSlot = enum {
    kernel,
    beta,
    incopy,
    itcopy,
    oncopy,
    otcopy,
    small_matrix_permit,
    small_kernel_nn,
    small_kernel_nt,
    small_kernel_nr,
    small_kernel_nc,
    small_kernel_tn,
    small_kernel_tt,
    small_kernel_tr,
    small_kernel_tc,
    small_kernel_rn,
    small_kernel_rt,
    small_kernel_rr,
    small_kernel_rc,
    small_kernel_cn,
    small_kernel_ct,
    small_kernel_cr,
    small_kernel_cc,
    small_kernel_b0_nn,
    small_kernel_b0_nt,
    small_kernel_b0_nr,
    small_kernel_b0_nc,
    small_kernel_b0_tn,
    small_kernel_b0_tt,
    small_kernel_b0_tr,
    small_kernel_b0_tc,
    small_kernel_b0_rn,
    small_kernel_b0_rt,
    small_kernel_b0_rr,
    small_kernel_b0_rc,
    small_kernel_b0_cn,
    small_kernel_b0_ct,
    small_kernel_b0_cr,
    small_kernel_b0_cc,

    /// The small-matrix kernel for the reference's `transa`, `transb` codes
    /// (`N` 0, `T` 1, `R` 2, `C` 3), the `beta = 0` one for `b0`.
    pub fn small(comptime b0: bool, comptime xa: u2, comptime xb: u2) GemmSlot {
        const first: GemmSlot = if (b0) .small_kernel_b0_nn else .small_kernel_nn;
        return @fromBackingInt(@intCast(@backingInt(first) + 4 * @as(usize, xa) + xb));
    }
};

/// The vtable field of `g`'s slot `s`: `sbgemm_kernel`.
pub fn gemmSlotName(comptime g: Gemm, comptime s: GemmSlot) []const u8 {
    return comptime g.prefix() ++ "_" ++ @tagName(s);
}

/// `gemmSlotName`'s `Slot`, for `blas.kernel`.
pub fn gemmSlotTag(comptime g: Gemm, comptime s: GemmSlot) util.Slot {
    return comptime @field(util.Slot, gemmSlotName(g, s));
}

/// `vt`'s slot `s` of `g`, by pointer; `vt` points to a `Vtable`.
pub inline fn gemmSlotPtr(vt: anytype, comptime g: Gemm, comptime s: GemmSlot) @TypeOf(&@field(vt, gemmSlotName(g, s))) {
    return &@field(vt, gemmSlotName(g, s));
}

/// The vtable field of half-precision `g`'s gemv (`sh`, `b`, `sb`), the
/// transposed one for `trans`: `sbgemv_t`.
pub fn halfGemvName(comptime g: Gemm, comptime trans: bool) []const u8 {
    return comptime @tagName(g) ++ "gemv_" ++ if (trans) "t" else "n";
}

/// `halfGemvName`'s `Slot`, for `blas.kernel`.
pub fn halfGemvSlot(comptime g: Gemm, comptime trans: bool) util.Slot {
    return comptime @field(util.Slot, halfGemvName(g, trans));
}

/// `vt`'s half-precision gemv field, by pointer; `vt` points to a `Vtable`.
pub inline fn halfGemvPtr(vt: anytype, comptime g: Gemm, comptime trans: bool) @TypeOf(&@field(vt, halfGemvName(g, trans))) {
    return &@field(vt, halfGemvName(g, trans));
}
