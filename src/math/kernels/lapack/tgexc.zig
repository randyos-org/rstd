//! `xTGEXC`: moving one diagonal block of a generalized Schur form,
//! Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dtgexc.f`, `ztgexc.f`,
//! OpenBLAS v0.3.34-317; modified BSD), over `tgex2.zig`'s swaps.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const schur = @import("schur.zig");
const tgex2 = @import("tgex2.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const isComplex = types.isComplex;
const Error = ops.Error;
const ReorderVectors = schur.ReorderVectors;

/// What `tgexc` returns: a slice too short or a leading dimension too small
/// (`q`, `z` checked only when updated), `RowOutOfRange` for an `ifst` or
/// `ilst` past the rows, `SwapRejected` for the reference's `INFO = 1`: two
/// adjacent blocks were too close to swap stably, `(A, B)` part reordered
/// and `ilst` the moving block's current row.
pub const TgexcError = interface.ShapeError || error{ RowOutOfRange, SwapRejected } || Allocator.Error || Io.Cancelable;

/// The pencil and its vectors as the moves take them.
pub fn Pencil(comptime T: type) type {
    return struct {
        n: usize,
        a: [*]T,
        lda: usize,
        b: [*]T,
        ldb: usize,
        q: ?[*]T,
        ldq: usize,
        z: ?[*]T,
        ldz: usize,

        fn qp(p: @This()) [*]T {
            return p.q orelse undefined;
        }

        fn zp(p: @This()) [*]T {
            return p.z orelse undefined;
        }
    };
}

/// `?TGEXC` (`dtgexc.f`, `ztgexc.f`): moves the diagonal block of the `n x
/// n` generalized Schur form `(A, B)` at row `ifst` (0-based) to row `ilst`
/// by orthogonal (unitary) equivalence, `(A, B) := Qᴴ (A, B) Z`, and as
/// `compq` and `compz` ask `Q := Q Qₘ`, `Z := Z Zₘ`. Real: `A` is
/// quasi-triangular with standardized 2 x 2 blocks; `ifst` moves to its
/// block's first row and `ilst` to the moved block's final first row, as the
/// reference's. Complex: both unchanged on success (the reference's
/// forward move reports one row short).
pub fn tgexc(allocator: Allocator, io: Io, comptime T: type, compq: ReorderVectors, compz: ReorderVectors, n: usize, a: []T, lda: usize, b: []T, ldb: usize, q: []T, ldq: usize, z: []T, ldz: usize, ifst: *usize, ilst: *usize) TgexcError!void {
    comptime ops.requireFloat(T, "tgexc");
    try util.holds(a.len, n, n, lda);
    try util.holds(b.len, n, n, ldb);
    if (compq == .update) try util.holds(q.len, n, n, ldq);
    if (compz == .update) try util.holds(z.len, n, n, ldz);
    if (n > 0 and (ifst.* >= n or ilst.* >= n)) return error.RowOutOfRange;
    if (n <= 1) return;
    const p: Pencil(T) = .{ .n = n, .a = a.ptr, .lda = lda, .b = b.ptr, .ldb = ldb, .q = if (compq == .update) q.ptr else null, .ldq = ldq, .z = if (compz == .update) z.ptr else null, .ldz = ldz };
    if (comptime isComplex(T)) {
        if (moveComplex(T, p, ifst.*, ilst.*)) |here| {
            ilst.* = here;
            return error.SwapRejected;
        }
        return;
    }
    // the reference's 1-based rows
    var f = ifst.* + 1;
    var l = ilst.* + 1;
    defer {
        ifst.* = f - 1;
        ilst.* = l - 1;
    }
    if (!try moveReal(allocator, io, T, p, &f, &l)) return error.SwapRejected;
}

/// `ztgexc.f` after its checks (0-based rows): null, or the moving block's
/// row when a swap failed.
pub fn moveComplex(comptime T: type, p: Pencil(T), ifst: usize, ilst: usize) ?usize {
    const swap = struct {
        fn f(pp: Pencil(T), j: usize) bool {
            return tgex2.tgex2Complex(T, pp.q != null, pp.z != null, pp.n, pp.a, pp.lda, pp.b, pp.ldb, pp.qp(), pp.ldq, pp.zp(), pp.ldz, j);
        }
    }.f;
    if (ifst < ilst) {
        var here = ifst;
        while (here < ilst) : (here += 1) if (!swap(p, here)) return here;
    } else {
        var here = ifst;
        while (here > ilst) : (here -= 1) if (!swap(p, here - 1)) return here;
    }
    return null;
}

/// `dtgexc.f` after its checks: `ifst` and `ilst` (1-based) as the
/// reference leaves them; false for `INFO = 1`.
pub fn moveReal(allocator: Allocator, io: Io, comptime F: type, p: Pencil(F), ifst: *usize, ilst: *usize) Error!bool {
    const n = p.n;
    // `A(i, j)` is nonzero, 1-based
    const nz = struct {
        fn f(pp: Pencil(F), i: usize, j: usize) bool {
            return pp.a[(i - 1) + (j - 1) * pp.lda] != 0;
        }
    }.f;
    // `DTGEX2` at 1-based `j1`
    const ex = struct {
        fn f(al: Allocator, i: Io, pp: Pencil(F), j1: usize, n1: usize, n2: usize) Error!bool {
            return tgex2.tgex2Real(al, i, F, pp.q != null, pp.z != null, pp.n, pp.a, pp.lda, pp.b, pp.ldb, pp.qp(), pp.ldq, pp.zp(), pp.ldz, j1 - 1, n1, n2);
        }
    }.f;
    if (ifst.* > 1 and nz(p, ifst.*, ifst.* - 1)) ifst.* -= 1;
    var nbf: usize = 1;
    if (ifst.* < n and nz(p, ifst.* + 1, ifst.*)) nbf = 2;
    if (ilst.* > 1 and nz(p, ilst.*, ilst.* - 1)) ilst.* -= 1;
    var nbl: usize = 1;
    if (ilst.* < n and nz(p, ilst.* + 1, ilst.*)) nbl = 2;
    if (ifst.* == ilst.*) return true;
    var here = ifst.*;
    if (ifst.* < ilst.*) {
        if (nbf == 2 and nbl == 1) ilst.* -= 1;
        if (nbf == 1 and nbl == 2) ilst.* += 1;
        while (true) {
            if (nbf == 1 or nbf == 2) {
                // the block is whole
                var nbnext: usize = 1;
                if (here + nbf + 1 <= n and nz(p, here + nbf + 1, here + nbf)) nbnext = 2;
                if (!try ex(allocator, io, p, here, nbf, nbnext)) {
                    ilst.* = here;
                    return false;
                }
                here += nbnext;
                if (nbf == 2 and !nz(p, here + 1, here)) nbf = 3;
            } else {
                // a 2 x 2 block split in two 1 x 1s: move each
                var nbnext: usize = 1;
                if (here + 3 <= n and nz(p, here + 3, here + 2)) nbnext = 2;
                if (!try ex(allocator, io, p, here + 1, 1, nbnext)) {
                    ilst.* = here;
                    return false;
                }
                if (nbnext == 1) {
                    if (!try ex(allocator, io, p, here, 1, 1)) {
                        ilst.* = here;
                        return false;
                    }
                    here += 1;
                } else {
                    if (!nz(p, here + 2, here + 1)) nbnext = 1;
                    if (nbnext == 2) {
                        if (!try ex(allocator, io, p, here, 1, nbnext)) {
                            ilst.* = here;
                            return false;
                        }
                        here += 2;
                    } else {
                        if (!try ex(allocator, io, p, here, 1, 1)) {
                            ilst.* = here;
                            return false;
                        }
                        here += 1;
                        if (!try ex(allocator, io, p, here, 1, 1)) {
                            ilst.* = here;
                            return false;
                        }
                        here += 1;
                    }
                }
            }
            if (!(here < ilst.*)) break;
        }
    } else {
        while (true) {
            if (nbf == 1 or nbf == 2) {
                var nbnext: usize = 1;
                if (here >= 3 and nz(p, here - 1, here - 2)) nbnext = 2;
                if (!try ex(allocator, io, p, here - nbnext, nbnext, nbf)) {
                    ilst.* = here;
                    return false;
                }
                here -= nbnext;
                if (nbf == 2 and !nz(p, here + 1, here)) nbf = 3;
            } else {
                var nbnext: usize = 1;
                if (here >= 3 and nz(p, here - 1, here - 2)) nbnext = 2;
                if (!try ex(allocator, io, p, here - nbnext, nbnext, 1)) {
                    ilst.* = here;
                    return false;
                }
                if (nbnext == 1) {
                    if (!try ex(allocator, io, p, here, nbnext, 1)) {
                        ilst.* = here;
                        return false;
                    }
                    here -= 1;
                } else {
                    if (!nz(p, here, here - 1)) nbnext = 1;
                    if (nbnext == 2) {
                        if (!try ex(allocator, io, p, here - 1, 2, 1)) {
                            ilst.* = here;
                            return false;
                        }
                        here -= 2;
                    } else {
                        if (!try ex(allocator, io, p, here, 1, 1)) {
                            ilst.* = here;
                            return false;
                        }
                        here -= 1;
                        if (!try ex(allocator, io, p, here, 1, 1)) {
                            ilst.* = here;
                            return false;
                        }
                        here -= 1;
                    }
                }
            }
            if (!(here > ilst.*)) break;
        }
    }
    ilst.* = here;
    return true;
}
