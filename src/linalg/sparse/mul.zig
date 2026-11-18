//! `A x` and `Aᵀ y` on CSC (src/linalg/design.md, "Sparse matrices", Threads),
//! and `A B` with `B` dense. `A x` scatters one column at a time; threaded,
//! its columns split into chunks of about equal nonzeros, each chunk's sum
//! in a partial vector over the rows it reaches, added to `y` in chunk order
//! (`threading.Partials`). `Aᵀ y` gathers each output whole, so its chunks
//! cannot change a bit. Chunks follow the size and pattern alone, so the bits
//! do too, at any worker count.

const std = @import("std");
const matrix = @import("../matrix.zig");
const types = @import("../blas/types.zig");
const threading = @import("../blas/threading.zig");
const ops = @import("../lapack/ops.zig");
const util = @import("util.zig");
const csc = @import("csc.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Trans = types.Trans;
const Chunk = threading.Chunk;

/// The routine family that sizes the work: a sparse `gemv`, `nnz` its work.
const family: threading.Family = .gemv;

pub fn mulVec(comptime T: type, a: csc.Csc(T), allocator: Allocator, io: Io, trans: Trans, alpha: T, x: []const T, beta: T, y: []T) csc.MulError!void {
    const tr = trans.transposed();
    if (x.len < @as(usize, if (tr) a.rows else a.cols) or y.len < @as(usize, if (tr) a.cols else a.rows)) return error.DimensionMismatch;
    switch (trans) {
        inline else => |t| {
            const g: Gemv(T, t) = .{ .a = a, .alpha = alpha, .x = x, .beta = beta, .y = y };
            return if (comptime t.transposed()) g.gather(io) else g.scatter(allocator, io);
        },
    }
}

/// `y := alpha op(A) x + beta y` for one `trans`, known at compile time.
fn Gemv(comptime T: type, comptime trans: Trans) type {
    return struct {
        a: csc.Csc(T),
        alpha: T,
        x: []const T,
        beta: T,
        y: []T,
        parts: ?*const threading.Partials(T) = null,

        const Self = @This();

        fn opA(v: T) T {
            return if (comptime trans.conjA()) ops.conj(T, v) else v;
        }

        fn opX(v: T) T {
            return if (comptime trans.conjX()) ops.conj(T, v) else v;
        }

        /// `y := beta y`, then each column's `alpha x_j` times it added in.
        fn scatter(g: Self, allocator: Allocator, io: Io) csc.MulError!void {
            const m: usize = g.a.rows;
            const y = g.y[0..m];
            if (ops.isZero(T, g.beta)) {
                @memset(y, ops.zero(T));
            } else {
                for (y) |*v| v.* = ops.mul(T, g.beta, v.*);
            }
            if (ops.isZero(T, g.alpha)) return;
            const nnz = g.a.nonzeros();
            if (threading.serial(family, T, nnz)) return g.columns(0, g.a.cols, y, 0);
            const p0 = threading.plan(family, T, nnz, g.a.cols, .even, 1);
            // each chunk at least `m` nonzeros, so its partial costs no more than its work
            const want = @min(p0.chunks.count, @max(nnz / @max(m, 1), 1));
            if (want == 1) return g.columns(0, g.a.cols, y, 0);
            const p: threading.Plan = .{ .chunks = split(T, g.a, want), .workers = p0.workers };
            const parts: threading.Partials(T) = try .init(allocator, &p, g, reach);
            defer parts.deinit(allocator);
            var h = g;
            h.parts = &parts;
            try threading.forEach(io, &p, @as(*const Self, &h), chunk);
            try parts.fold(io, y.ptr, 1);
        }

        /// The rows chunk `c`'s columns reach: from its least first row to
        /// its greatest last, rows being ascending.
        fn reach(g: Self, c: Chunk) threading.Rows {
            var lo: usize = std.math.maxInt(usize);
            var hi: usize = 0;
            for (c.start..c.end) |j| {
                const s = g.a.col_start[j];
                const e = g.a.col_start[j + 1];
                if (s == e) continue;
                lo = @min(lo, g.a.row[s]);
                hi = @max(hi, @as(usize, g.a.row[e - 1]) + 1);
            }
            return if (hi == 0) .{ .lo = 0, .hi = 0 } else .{ .lo = lo, .hi = hi };
        }

        fn chunk(g: *const Self, c: Chunk) void {
            const parts = g.parts.?;
            g.columns(c.start, c.end, parts.take(c.index), parts.rows[c.index].lo);
        }

        /// Columns `j0..j1` added into `out`, whose element 0 is row `lo`.
        fn columns(g: Self, j0: usize, j1: usize, out: []T, lo: usize) void {
            for (j0..j1) |j| {
                const t = ops.mul(T, g.alpha, opX(g.x[j]));
                for (g.a.col_start[j]..g.a.col_start[j + 1]) |k| {
                    const r = g.a.row[k] - lo;
                    out[r] = ops.add(T, out[r], ops.mul(T, opA(g.a.value[k]), t));
                }
            }
        }

        /// Each `y_j := alpha op(A)(:, j) x + beta y_j`, its dot whole.
        fn gather(g: Self, io: Io) csc.MulError!void {
            const nnz = g.a.nonzeros();
            if (threading.serial(family, T, nnz)) return g.dots(.{ .index = 0, .start = 0, .end = g.a.cols });
            const p0 = threading.plan(family, T, nnz, g.a.cols, .even, 1);
            const p: threading.Plan = .{ .chunks = split(T, g.a, p0.chunks.count), .workers = p0.workers };
            return threading.forEach(io, &p, g, dots);
        }

        fn dots(g: Self, c: Chunk) void {
            for (c.start..c.end) |j| {
                var s = ops.zero(T);
                if (!ops.isZero(T, g.alpha)) {
                    for (g.a.col_start[j]..g.a.col_start[j + 1]) |k| s = ops.add(T, s, ops.mul(T, opA(g.a.value[k]), opX(g.x[g.a.row[k]])));
                }
                const as = ops.mul(T, g.alpha, s);
                g.y[j] = if (ops.isZero(T, g.beta)) as else ops.add(T, as, ops.mul(T, g.beta, g.y[j]));
            }
        }
    };
}

/// `a`'s columns in up to `want` chunks of about equal nonzeros, empty ones
/// dropped; one empty chunk without columns.
fn split(comptime T: type, a: csc.Csc(T), want: usize) threading.Chunks {
    const c = std.math.clamp(want, 1, threading.chunks_max);
    const nnz: u64 = a.nonzeros();
    var r: threading.Chunks = .{ .count = 0, .bounds = undefined };
    r.bounds[0] = 0;
    for (1..c + 1) |k| {
        const b = if (k == c) a.cols else util.lowerBound(a.col_start, nnz * k / c);
        if (b > r.bounds[r.count]) {
            r.count += 1;
            r.bounds[r.count] = b;
        }
    }
    if (r.count == 0) {
        r.count = 1;
        r.bounds[1] = 0;
    }
    return r;
}

/// `C := A B`, each column of `C` scattered whole, so the chunks of `B`'s
/// columns cannot change a bit.
pub fn mulDense(comptime T: type, a: csc.Csc(T), allocator: Allocator, io: Io, b: matrix.ConstView(T)) csc.MulError!matrix.Mat(T) {
    if (b.rows != a.cols) return error.DimensionMismatch;
    const c: matrix.Mat(T) = try .zeros(allocator, a.rows, b.cols);
    errdefer c.deinit(allocator);
    const D = struct {
        a: csc.Csc(T),
        b: matrix.ConstView(T),
        c: matrix.Mat(T),

        fn cols(d: @This(), ch: Chunk) void {
            const m: usize = d.a.rows;
            for (ch.start..ch.end) |q| {
                const out = d.c.data[q * m ..][0..m];
                for (d.b.col(q), 0..) |bj, j| {
                    for (d.a.col_start[j]..d.a.col_start[j + 1]) |k| out[d.a.row[k]] = ops.add(T, out[d.a.row[k]], ops.mul(T, d.a.value[k], bj));
                }
            }
        }
    };
    const d: D = .{ .a = a, .b = b, .c = c };
    const work = threading.product(&.{ a.nonzeros(), b.cols });
    if (threading.serial(family, T, work)) {
        d.cols(.{ .index = 0, .start = 0, .end = b.cols });
    } else {
        const p = threading.plan(family, T, work, b.cols, .even, 1);
        try threading.forEach(io, &p, d, D.cols);
    }
    return c;
}
