//! The portable `omatcopy` and `imatcopy` kernels, column-major: `B :=
//! alpha*op(A)` out of place, and in place where `op` keeps the shape.
//! `kernel/arm/omatcopy_cn.c`, `omatcopy_ct.c`, `zomatcopy_cn.c`,
//! `zomatcopy_cnc.c`, `zomatcopy_ct.c`, `zomatcopy_ctc.c`, and
//! `kernel/generic/imatcopy_cn.c`, `imatcopy_ct.c`, `zimatcopy_cn.c`,
//! `zimatcopy_cnc.c`, `zimatcopy_ct.c`, `zimatcopy_ctc.c`; the row-major
//! `_r*` ones are left out.
//!
//! Every kernel scales an element alike (`scale`): a zero `alpha` stores
//! zeros without reading `A`, and `alpha` 1 copies (or conjugates) its bits;
//! else the products of the C. Upstream's real kernels do both; the complex
//! ones multiply always, so `alpha` 1 there turns `0 * Inf` into NaN in the
//! other part, and `zimatcopy`'s second copy, by 1, does so even for another
//! `alpha` (`issues.md`). One rule for all keeps the in-place and buffered
//! `imatcopy` and `omatcopy` bit for bit the same.

const std = @import("std");
const types = @import("../types.zig");
const util = @import("../util.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const check = @import("../_test/check.zig");
const bitwise = @import("../_test/bitwise.zig");
const pattern = @import("../_test/pattern.zig");

const BlasLong = types.BlasLong;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;

/// `alpha*x`, `x` conjugated first where `conj`; `T`'s element as the kernels
/// see it, two floats for a complex one. Complex by `util.cmul`, so every
/// kernel's bits match on the x87 too.
fn Scale(comptime T: type) type {
    const complex = types.isComplex(T);
    const F = types.Float(T);
    return struct {
        const E = if (complex) [2]F else F;
        const Alpha = if (complex) [2]F else F;

        fn of(alpha: Alpha, comptime conj: bool, x: E) E {
            if (comptime !complex) return if (alpha == 0) 0 else if (alpha == 1) x else alpha * x;
            const ar = alpha[0];
            const ai = alpha[1];
            if (ar == 0 and ai == 0) return .{ 0, 0 };
            if (ar == 1 and ai == 0) return if (conj) .{ x[0], -x[1] } else x;
            return util.cmul(F, ar, ai, x[0], if (conj) -x[1] else x[1]);
        }

        /// `B := alpha*op(A)`, `A` `rows x cols`; `op` transposes where
        /// `trans`, so `B` is `cols x rows` then.
        fn outOf(alpha: Alpha, comptime trans: bool, comptime conj: bool, rows: BlasLong, cols: BlasLong, a: [*]const E, lda: BlasLong, b: [*]E, ldb: BlasLong) void {
            if (rows <= 0 or cols <= 0) return;
            const ia: usize = @intCast(lda);
            const ib: usize = @intCast(ldb);
            for (0..@intCast(cols)) |j| for (0..@intCast(rows)) |i| {
                const v = of(alpha, conj, a[i + j * ia]);
                if (trans) b[j + i * ib] = v else b[i + j * ib] = v;
            };
        }

        /// `A := alpha*op(A)` in place; transposed only for a square `A`.
        fn inPlace(alpha: Alpha, comptime trans: bool, comptime conj: bool, rows: BlasLong, cols: BlasLong, a: [*]E, lda: BlasLong) void {
            if (rows <= 0 or cols <= 0) return;
            const ia: usize = @intCast(lda);
            if (!trans) {
                for (0..@intCast(cols)) |j| for (0..@intCast(rows)) |i| {
                    a[i + j * ia] = of(alpha, conj, a[i + j * ia]);
                };
                return;
            }
            const n: usize = @intCast(@min(rows, cols));
            for (0..n) |j| {
                a[j + j * ia] = of(alpha, conj, a[j + j * ia]);
                for (j + 1..n) |i| {
                    const lo = a[i + j * ia];
                    a[i + j * ia] = of(alpha, conj, a[j + i * ia]);
                    a[j + i * ia] = of(alpha, conj, lo);
                }
            }
        }
    };
}

/// The kernels over `T`, by slot less its letter: `omatcopy_k_cn` and the
/// rest, `cnc` and `ctc` conjugating (complex only).
fn Kernels(comptime T: type) type {
    const F = types.Float(T);
    const cc = types.kernelCc(F);
    const S = Scale(T);
    return if (types.isComplex(T)) struct {
        fn O(comptime trans: bool, comptime conj: bool) type {
            return struct {
                fn k(rows: BlasLong, cols: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong, b: [*]F, ldb: BlasLong) callconv(cc) i32 {
                    S.outOf(.{ alpha_r, alpha_i }, trans, conj, rows, cols, @ptrCast(a), lda, @ptrCast(b), ldb);
                    return 0;
                }
            };
        }
        fn I(comptime trans: bool, comptime conj: bool) type {
            return struct {
                fn k(rows: BlasLong, cols: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong) callconv(cc) i32 {
                    S.inPlace(.{ alpha_r, alpha_i }, trans, conj, rows, cols, @ptrCast(a), lda);
                    return 0;
                }
            };
        }
        pub const omatcopy_k_cn = &O(false, false).k;
        pub const omatcopy_k_ct = &O(true, false).k;
        pub const omatcopy_k_cnc = &O(false, true).k;
        pub const omatcopy_k_ctc = &O(true, true).k;
        pub const imatcopy_k_cn = &I(false, false).k;
        pub const imatcopy_k_ct = &I(true, false).k;
        pub const imatcopy_k_cnc = &I(false, true).k;
        pub const imatcopy_k_ctc = &I(true, true).k;
    } else struct {
        fn O(comptime trans: bool) type {
            return struct {
                fn k(rows: BlasLong, cols: BlasLong, alpha: F, a: [*]F, lda: BlasLong, b: [*]F, ldb: BlasLong) callconv(cc) i32 {
                    S.outOf(alpha, trans, false, rows, cols, a, lda, b, ldb);
                    return 0;
                }
            };
        }
        fn I(comptime trans: bool) type {
            return struct {
                fn k(rows: BlasLong, cols: BlasLong, alpha: F, a: [*]F, lda: BlasLong) callconv(cc) i32 {
                    S.inPlace(alpha, trans, false, rows, cols, a, lda);
                    return 0;
                }
            };
        }
        pub const omatcopy_k_cn = &O(false).k;
        pub const omatcopy_k_ct = &O(true).k;
        pub const imatcopy_k_cn = &I(false).k;
        pub const imatcopy_k_ct = &I(true).k;
    };
}

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub const Portable = Kernels;

/// The precisions with slots; upstream has no `q` or `x` ones.
const precisions = .{ .{ "s", f32 }, .{ "d", f64 }, .{ "c", Complex32 }, .{ "z", Complex64 } };

/// Points every `omatcopy` and `imatcopy` slot at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    inline for (precisions) |p| {
        const K = Kernels(p[1]);
        inline for (comptime std.meta.declarations(K)) |d| {
            if (comptime @hasField(Vtable, p[0] ++ d)) @field(vt, p[0] ++ d) = @field(K, d);
        }
    }
}

/// A kernel by its suffix: `t` transposes, a trailing `c` conjugates (complex only).
pub const Copy = enum {
    cn,
    ct,
    cnc,
    ctc,

    /// The `omatcopy` kernel's op.
    pub fn outOp(c: Copy) portable.Op {
        return switch (c) {
            .cn => .omatcopy_k_cn,
            .ct => .omatcopy_k_ct,
            .cnc => .omatcopy_k_cnc,
            .ctc => .omatcopy_k_ctc,
        };
    }

    /// The `imatcopy` kernel's op.
    pub fn inOp(c: Copy) portable.Op {
        return switch (c) {
            .cn => .imatcopy_k_cn,
            .ct => .imatcopy_k_ct,
            .cnc => .imatcopy_k_cnc,
            .ctc => .imatcopy_k_ctc,
        };
    }
};

// ---- tests -----------------------------------------------------------------

/// `alpha*op(A)(i, j)` by the definition, in `T`.
fn want(comptime T: type, op: Copy, alpha: T, a: []const T, lda: usize, i: usize, j: usize) T {
    const t = op == .ct or op == .ctc;
    var x = if (t) a[j + i * lda] else a[i + j * lda];
    if (comptime types.isComplex(T)) {
        if (op == .cnc or op == .ctc) x = x.conjugate();
        if (alpha.re == 0 and alpha.im == 0) return .init(0, 0);
        if (alpha.re == 1 and alpha.im == 0) return x;
        const p = util.cmul(types.Float(T), alpha.re, alpha.im, x.re, x.im);
        return .init(p[0], p[1]);
    }
    return if (alpha == 0) 0 else alpha * x;
}

/// Every kernel of `vt` over `T` against `want`, bit for bit, out of place and
/// in place, with special values in `A` and in the padding past its rows.
fn checkType(comptime T: type, comptime pre: []const u8, vt: *const Vtable) !void {
    const complex = comptime types.isComplex(T);
    const F = types.Float(T);
    const al = std.testing.allocator;
    const ops: []const Copy = comptime if (complex) &.{ .cn, .ct, .cnc, .ctc } else &.{ .cn, .ct };
    const alphas = if (complex) [_]T{ T.init(1.5, -0.75), T.init(0, 0), T.init(1, 0), T.init(0, 1) } else [_]T{ 1.5, 0, 1, -1 };
    const sizes = [_][2]usize{ .{ 1, 1 }, .{ 3, 1 }, .{ 1, 4 }, .{ 5, 3 }, .{ 7, 7 }, .{ 17, 9 } };
    inline for (ops) |op| for (alphas) |alpha| for (sizes) |s| {
        const rows = s[0];
        const cols = s[1];
        const t = op == .ct or op == .ctc;
        const lda = rows + 2;
        const ldb = (if (t) cols else rows) + 1;
        const a = try al.alloc(T, lda * cols);
        defer al.free(a);
        const b = try al.alloc(T, ldb * (if (t) rows else cols));
        defer al.free(b);
        check.fill(T, a, rows * 13 + cols);
        for (a, 0..) |*v, i| v.* = pattern.shape(T, v.*, i, .mixed);
        const sentinel: T = if (complex) T.init(7, 7) else 7;
        @memset(b, sentinel);
        const ar: [*]F = @ptrCast(a.ptr);
        const br: [*]F = @ptrCast(b.ptr);
        const k = portable.slotPtr(vt, .from(pre), op.outOp()).*.?;
        if (complex) {
            _ = k(@intCast(rows), @intCast(cols), alpha.re, alpha.im, ar, @intCast(lda), br, @intCast(ldb));
        } else {
            _ = k(@intCast(rows), @intCast(cols), alpha, ar, @intCast(lda), br, @intCast(ldb));
        }
        const out_rows = if (t) cols else rows;
        const out_cols = if (t) rows else cols;
        for (0..out_cols) |j| for (0..ldb) |i| {
            const got = b[i + j * ldb];
            const w = if (i < out_rows) want(T, op, alpha, a, lda, i, j) else sentinel;
            if (!same(T, w, got)) {
                check.say("{s}omatcopy_k_{s} {d}x{d} alpha {any} at ({d},{d}): ", .{ pre, @tagName(op), rows, cols, alpha, i, j });
                return error.TestUnexpectedResult;
            }
        };
        // in place: the shape kept, or square; the result is the same
        if (t and rows != cols) continue;
        const ik = portable.slotPtr(vt, .from(pre), op.inOp()).*.?;
        if (complex) {
            _ = ik(@intCast(rows), @intCast(cols), alpha.re, alpha.im, ar, @intCast(lda));
        } else {
            _ = ik(@intCast(rows), @intCast(cols), alpha, ar, @intCast(lda));
        }
        for (0..out_cols) |j| for (0..out_rows) |i| {
            if (!bitwise.same(T, b[i + j * ldb], a[i + j * lda], true)) {
                check.say("{s}imatcopy_k_{s} {d}x{d} alpha {any} at ({d},{d}): ", .{ pre, @tagName(op), rows, cols, alpha, i, j });
                return error.TestUnexpectedResult;
            }
        };
    };
}

/// `bitwise.same`, part by part, but any NaN for a NaN: which operand's
/// payload a product keeps varies with the target.
fn same(comptime T: type, w: T, g: T) bool {
    if (comptime types.isComplex(T)) return same(types.Float(T), w.re, g.re) and same(types.Float(T), w.im, g.im);
    return bitwise.same(T, w, g, true) or (std.math.isNan(w) and std.math.isNan(g));
}

test "portable omatcopy and imatcopy kernels follow the definition, and agree with each other bit for bit" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    inline for (precisions) |p| try checkType(p[1], p[0], &vt);
}
