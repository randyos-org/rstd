//! `trmm` and `trsm`: OpenBLAS's `driver/level3/trmm_L.c`, `trmm_R.c`,
//! `trsm_L.c` and `trsm_R.c` on `driver.zig`'s blocking, packing space and
//! kernel calls; more than one worker splits `B` into serial parts.

const std = @import("std");
const types = @import("../types.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const util = @import("../util.zig");
const threading = @import("../threading.zig");
const driver = @import("driver.zig");
const tmm = @import("triangular/trmm.zig");
const tsm = @import("triangular/trsm.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const isComplex = types.isComplex;
const Op = portable.Op;
const Side = driver.Side;
const Blocking = driver.Blocking;
const Conj = driver.Conj;
const Error = driver.Error;
const Work = driver.Work;
const scale = driver.scale;
const unlessOne = driver.unlessOne;
const placed = driver.placed;
const stack_bytes = driver.stack_bytes;
const raw = util.raw;
const isZero = util.isZero;
const trmmL = tmm.trmmL;
const trmmR = tmm.trmmR;
const trsmL = tsm.trsmL;
const trsmR = tsm.trsmR;

/// A copy of the triangle (`trsm_iunucopy` and the like); `offset` places
/// the block's lines against its depth.
fn triCopy(comptime T: type, vt: *const Vtable, comptime op: Op, m: usize, n: usize, a: [*]T, lda: usize, offset: isize, dst: [*]T) void {
    if (comptime isComplex(T)) {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), raw(T, a), @intCast(lda), offset, raw(T, dst));
    } else {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), a, @intCast(lda), offset, dst);
    }
}

/// A `trmm` or `trsm` kernel: `alpha` for `trmm`; `trsm`'s, a dummy
/// upstream, is `unitTag`.
fn triKernel(comptime T: type, vt: *const Vtable, comptime op: Op, m: usize, n: usize, k: usize, alpha: T, sa: [*]T, sb: [*]T, c: [*]T, ldc: usize, offset: isize) void {
    if (comptime isComplex(T)) {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), @intCast(k), alpha.re, alpha.im, raw(T, sa), raw(T, sb), raw(T, c), @intCast(ldc), offset);
    } else {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), @intCast(k), alpha, sa, sb, c, @intCast(ldc), offset);
    }
}

pub const fromReal = util.fromReal;

/// The `trsm` kernels' `alpha`: 1 for a unit diagonal, which they skip
/// dividing by, else -1 (upstream's `dm1`).
pub fn unitTag(comptime T: type, unit: bool) T {
    return fromReal(T, if (unit) 1 else -1);
}

/// The options of a `trmm` or `trsm`: values, as `Desc`'s are.
pub const Tri = struct {
    left: bool,
    upper: bool,
    trans: bool,
    conj: bool,
    unit: bool,

    /// The gemm kernel's conjugation: the triangle's side.
    pub fn gemmConj(t: Tri) Conj {
        return if (!t.conj) .n else if (t.left) .l else .r;
    }
};

/// Which triangle copy a driver packs with: its triangle and whether it reads
/// across rows (`t`) or down columns (`n`).
pub const TriCopy = struct { upper: bool, across: bool };

/// A triangular routine, by its kernels' prefix.
const Routine = enum { trmm, trsm };

/// `<r>_<io><u|l><n|t><u|n>copy`: the triangle, across rows (`t`) or down
/// columns, and a unit diagonal.
fn copyOp(comptime r: Routine, comptime io: Side, comptime upper: bool, comptime across: bool, comptime unit: bool) Op {
    // by upper, across, unit
    const table: [2][2][2]Op = switch (r) {
        .trmm => switch (io) {
            .i => .{ .{ .{ .trmm_ilnncopy, .trmm_ilnucopy }, .{ .trmm_iltncopy, .trmm_iltucopy } }, .{ .{ .trmm_iunncopy, .trmm_iunucopy }, .{ .trmm_iutncopy, .trmm_iutucopy } } },
            .o => .{ .{ .{ .trmm_olnncopy, .trmm_olnucopy }, .{ .trmm_oltncopy, .trmm_oltucopy } }, .{ .{ .trmm_ounncopy, .trmm_ounucopy }, .{ .trmm_outncopy, .trmm_outucopy } } },
        },
        .trsm => switch (io) {
            .i => .{ .{ .{ .trsm_ilnncopy, .trsm_ilnucopy }, .{ .trsm_iltncopy, .trsm_iltucopy } }, .{ .{ .trsm_iunncopy, .trsm_iunucopy }, .{ .trsm_iutncopy, .trsm_iutucopy } } },
            .o => .{ .{ .{ .trsm_olnncopy, .trsm_olnucopy }, .{ .trsm_oltncopy, .trsm_oltucopy } }, .{ .{ .trsm_ounncopy, .trsm_ounucopy }, .{ .trsm_outncopy, .trsm_outucopy } } },
        },
    };
    return table[@intFromBool(upper)][@intFromBool(across)][@intFromBool(unit)];
}

/// `<r>_kernel_<L|R><N|T>`, or `R` and `C` for `N` and `T` conjugated.
fn kernelOp(comptime r: Routine, comptime left: bool, comptime nk: bool, comptime conj: bool) Op {
    // by left, conj, nk
    const table: [2][2][2]Op = switch (r) {
        .trmm => .{ .{ .{ .trmm_kernel_RT, .trmm_kernel_RN }, .{ .trmm_kernel_RC, .trmm_kernel_RR } }, .{ .{ .trmm_kernel_LT, .trmm_kernel_LN }, .{ .trmm_kernel_LC, .trmm_kernel_LR } } },
        .trsm => .{ .{ .{ .trsm_kernel_RT, .trsm_kernel_RN }, .{ .trsm_kernel_RC, .trsm_kernel_RR } }, .{ .{ .trsm_kernel_LT, .trsm_kernel_LN }, .{ .trsm_kernel_LC, .trsm_kernel_LR } } },
    };
    return table[@intFromBool(left)][@intFromBool(conj)][@intFromBool(nk)];
}

/// `copyOp`'s copy (`trmm`'s by `posX`, `posY`; `trsm`'s by `offset`), as
/// `c` and `unit` say.
pub fn triPack(comptime T: type, vt: *const Vtable, comptime r: Routine, comptime io: Side, c: TriCopy, unit: bool, m: usize, n: usize, a: [*]T, lda: usize, x: usize, y: usize, offset: isize, dst: [*]T) void {
    switch (c.upper) {
        inline else => |u| switch (c.across) {
            inline else => |t| switch (unit) {
                inline else => |un| {
                    const op = comptime copyOp(r, io, u, t, un);
                    if (comptime r == .trmm) return placed(T, vt, op, m, n, a, lda, x, y, dst);
                    return triCopy(T, vt, op, m, n, a, lda, offset, dst);
                },
            },
        },
    }
}

/// `kernelOp`'s kernel for `o`'s side and conjugation.
pub fn triRun(comptime T: type, vt: *const Vtable, comptime r: Routine, o: Tri, nk: bool, m: usize, n: usize, k: usize, alpha: T, sa: [*]T, sb: [*]T, c: [*]T, ldc: usize, offset: isize) void {
    switch (o.left) {
        inline else => |l| switch (nk) {
            inline else => |nn| switch (o.conj) {
                inline else => |cj| if (comptime cj and !isComplex(T)) unreachable else {
                    return triKernel(T, vt, comptime kernelOp(r, l, nn, cj), m, n, k, alpha, sa, sb, c, ldc, offset);
                },
            },
        },
    }
}

/// `op(A)` from row `i`, column `l`, as the copies take it.
pub fn opAt(comptime T: type, trans: bool, p: [*]T, ld: usize, i: usize, l: usize) [*]T {
    return if (trans) p + l + i * ld else p + i + l * ld;
}

/// The packing space of `trmm` and `trsm`, as `Work(T).lens` gives it: a row
/// block of `P` by a slice of `Q`, and that slice by `R` columns (the
/// reference's kernels take all of it).
pub fn triLens(blk: Blocking, m: usize, n: usize, ka: usize) [3]usize {
    const depth = @min(ka, blk.q);
    return .{ @min(m, blk.p) * depth, depth * @min(n, blk.r), 0 };
}

/// `B := alpha*op(A)*B` (`left`) or `alpha*B*op(A)`, `A` triangular
/// (`upper` or lower, `unit` diagonal), `op` transposing (`trans`) and
/// conjugating (`conj`, complex) as asked; `B` is `m x n`. `B` is scaled by
/// `alpha` first, never read for a zero one (`trmm_L.c`, `trmm_R.c`). Up to
/// `want` workers run it.
pub fn trmm(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, left: bool, upper: bool, trans: bool, conj: bool, unit: bool, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]T, ldb: usize, want: u32) Error!void {
    return triangular(allocator, io, T, vt, false, left, upper, trans, conj, unit, m, n, alpha, a, lda, b, ldb, want);
}

/// `B := alpha*inv(op(A))*B` (`left`) or `alpha*B*inv(op(A))`, as `trmm`
/// takes its options (`trsm_L.c`, `trsm_R.c`). Up to `want` workers run it.
pub fn trsm(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, left: bool, upper: bool, trans: bool, conj: bool, unit: bool, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]T, ldb: usize, want: u32) Error!void {
    return triangular(allocator, io, T, vt, true, left, upper, trans, conj, unit, m, n, alpha, a, lda, b, ldb, want);
}

fn triangular(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, comptime solve: bool, left: bool, upper: bool, trans: bool, conj: bool, unit: bool, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]T, ldb: usize, want: u32) Error!void {
    // the reference passes `alpha` as `beta`; the drivers fold it in
    if (isZero(T, alpha)) return scale(T, vt, m, n, alpha, b, ldb);
    const s = unlessOne(T, alpha);
    const blk = Blocking.of(T, vt);
    const ap: [*]T = @constCast(a);
    const opt: Tri = .{ .left = left, .upper = upper, .trans = trans, .conj = isComplex(T) and conj, .unit = unit };
    if (want > 1) return parallel(allocator, io, T, vt, solve, opt, blk, m, n, ap, lda, s, b, ldb, want);
    var stack: [stack_bytes]u8 align(Work(T).alignment) = undefined;
    const lens = triLens(blk, m, n, if (left) m else n);
    const sc = try util.Scratch(T, Work(T).alignment, Work(T).alignment).init(allocator, &stack, Work(T).total(lens));
    defer sc.deinit(allocator);
    triDrive(T, vt, solve, opt, blk, m, n, ap, lda, s, b, ldb, Work(T).at(sc.ptr, lens));
}

/// The `trmm` (or with `solve`, `trsm`) loops on `work`.
pub fn triDrive(comptime T: type, vt: *const Vtable, comptime solve: bool, o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, s: ?T, b: [*]T, ldb: usize, work: Work(T)) void {
    if (solve) {
        if (o.left) trsmL(T, vt, o, blk, m, n, a, lda, s, b, ldb, work) else trsmR(T, vt, o, blk, m, n, a, lda, s, b, ldb, work);
    } else {
        if (o.left) trmmL(T, vt, o, blk, m, n, a, lda, s, b, ldb, work) else trmmR(T, vt, o, blk, m, n, a, lda, s, b, ldb, work);
    }
}

/// `B := s*B` on `m x n` at `b`, unless `s` is null: the triangular
/// drivers' `alpha`, where each part of `B` is first read.
pub fn scaleIf(comptime T: type, vt: *const Vtable, s: ?T, m: usize, n: usize, b: [*]T, ldb: usize) void {
    if (s) |v| scale(T, vt, m, n, v, b, ldb);
}

// ---- threaded --------------------------------------------------------------------

/// `trmm` or, with `solve`, `trsm` on up to `want` workers: column
/// panels of `B` on the left, row panels on the right, as OpenBLAS's
/// `gemm_thread_n` and `gemm_thread_m` divide it and BLIS splits its `jc` and
/// `ic` loops. Each part is an independent serial problem on packing space
/// of its own, so they run through `threading.forEach`; parts start on
/// `unroll_n` columns or `unroll_m` rows, so tiles are the serial ones.
fn parallel(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, comptime solve: bool, o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, s: ?T, b: [*]T, ldb: usize, want: u32) Error!void {
    const W = Work(T);
    const Parts = struct {
        vt: *const Vtable,
        o: Tri,
        blk: Blocking,
        m: usize,
        n: usize,
        a: [*]T,
        lda: usize,
        s: ?T,
        b: [*]T,
        ldb: usize,

        /// Part `c`'s `m`, `n` and `B`.
        fn of(p: *const @This(), c: threading.Chunk) struct { usize, usize, [*]T } {
            return if (p.o.left) .{ p.m, c.end - c.start, p.b + c.start * p.ldb } else .{ c.end - c.start, p.n, p.b + c.start };
        }

        fn lens(p: *const @This(), c: threading.Chunk) [3]usize {
            const d = p.of(c);
            return triLens(p.blk, d[0], d[1], if (p.o.left) p.m else p.n);
        }

        fn len(p: *const @This(), c: threading.Chunk) usize {
            return W.total(p.lens(c));
        }

        fn one(p: *const @This(), c: threading.Chunk, buf: []T) void {
            const d = p.of(c);
            triDrive(T, p.vt, solve, p.o, p.blk, d[0], d[1], p.a, p.lda, p.s, d[2], p.ldb, W.at(buf.ptr, p.lens(c)));
        }
    };
    const plan: threading.Plan = .{
        .chunks = .init(if (o.left) n else m, want, .even, if (o.left) blk.un else blk.um),
        .workers = want,
    };
    const parts: Parts = .{ .vt = vt, .o = o, .blk = blk, .m = m, .n = n, .a = a, .lda = lda, .s = s, .b = b, .ldb = ldb };
    try threading.forEachPart(allocator, io, &plan, T, W.alignment, &parts, Parts.len, Parts.one);
}
