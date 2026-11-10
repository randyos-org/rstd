//! Vtable slot names by enum: the ops by precision, the gemms, their
//! blocking parameters and kernel slots, and the half-precision gemvs. A leaf,
//! as the kernel objects compile it (`ext/bf16/kernel.zig`, `vector.zig`);
//! `portable.zig` re-exports it.

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
    geadd_k,
    omatcopy_k_cn,
    omatcopy_k_ct,
    omatcopy_k_cnc,
    omatcopy_k_ctc,
    imatcopy_k_cn,
    imatcopy_k_ct,
    imatcopy_k_cnc,
    imatcopy_k_ctc,
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
    laswp_ncopy,
    gemm3m_kernel,
    gemm3m_incopyb,
    gemm3m_incopyr,
    gemm3m_incopyi,
    gemm3m_itcopyb,
    gemm3m_itcopyr,
    gemm3m_itcopyi,
    gemm3m_oncopyb,
    gemm3m_oncopyr,
    gemm3m_oncopyi,
    gemm3m_otcopyb,
    gemm3m_otcopyr,
    gemm3m_otcopyi,
    symm3m_iucopyb,
    symm3m_iucopyr,
    symm3m_iucopyi,
    symm3m_ilcopyb,
    symm3m_ilcopyr,
    symm3m_ilcopyi,
    symm3m_oucopyb,
    symm3m_oucopyr,
    symm3m_oucopyi,
    symm3m_olcopyb,
    symm3m_olcopyr,
    symm3m_olcopyi,
    hemm3m_iucopyb,
    hemm3m_iucopyr,
    hemm3m_iucopyi,
    hemm3m_ilcopyb,
    hemm3m_ilcopyr,
    hemm3m_ilcopyi,
    hemm3m_oucopyb,
    hemm3m_oucopyr,
    hemm3m_oucopyi,
    hemm3m_olcopyb,
    hemm3m_olcopyr,
    hemm3m_olcopyi,
};

/// A slot's precision, by its letter: `s`, `d` and `q` real, `c`, `z` and
/// `x` their complex. `q` and `x` are `XDouble`'s, which may be `f64`.
pub const Prec = enum {
    s,
    d,
    q,
    c,
    z,
    x,

    /// The precision a letter (`"s"`) names. `@field`, not `stringToEnum`,
    /// whose comptime map costs every caller's branch quota.
    pub fn from(comptime l: []const u8) Prec {
        return @field(Prec, l);
    }

    /// The real precision of a complex one, a real one itself.
    pub fn real(p: Prec) Prec {
        return switch (p) {
            .s, .c => .s,
            .d, .z => .d,
            .q, .x => .q,
        };
    }

    pub fn isComplex(p: Prec) bool {
        return p != p.real();
    }
};

/// The vtable field `op` fills at `prec`: `.s` and `.axpy_k` give `saxpy_k`,
/// `.c` and `.rot_k` `csrot_k`.
pub fn slotName(comptime prec: Prec, comptime op: Op) []const u8 {
    return comptime blk: {
        const l = @tagName(prec);
        const tag = @tagName(op);
        break :blk switch (op) {
            .iamax_k, .iamin_k, .imax_k, .imin_k => "i" ++ l ++ tag[1..],
            .rot_k => if (prec.isComplex()) l ++ @tagName(prec.real()) ++ tag else l ++ tag,
            else => l ++ tag,
        };
    };
}

/// The `Slot` of `op` at `prec`, for `blas.kernel`.
pub fn slotTag(comptime prec: Prec, comptime op: Op) util.Slot {
    return comptime @field(util.Slot, slotName(prec, op));
}

/// `vt`'s field for `op` at `prec`, by pointer; `vt` points to a `Vtable`.
pub inline fn slotPtr(vt: anytype, comptime prec: Prec, comptime op: Op) @TypeOf(&@field(vt, slotName(prec, op))) {
    return &@field(vt, slotName(prec, op));
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
