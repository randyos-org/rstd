//! The 3M drivers: `driver/level3/gemm3m_level3.c`, and through it
//! `symm3m_k.c` and `hemm3m_k.c`, on `driver.zig`'s loops: three real
//! products a slice where `gemm` takes four, normwise stable only (`gemm/notes.md`).

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../types.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const util = @import("../util.zig");
const threading = @import("../threading.zig");
const test_io = @import("../_test/io.zig");
const driver = @import("driver.zig");
const Part = @import("gemm3m/kernel.zig").Part;

const Allocator = std.mem.Allocator;
const Io = std.Io;
const BlasLong = types.BlasLong;
const Float = types.Float;
const Blocking = driver.Blocking;
const Matrix = driver.Matrix;
const Pack = driver.Pack;
const Op = portable.Op;

pub const Error = driver.Error;

/// `GEMM3M_P`, `_Q`, `_R` and the 3M unrolls from the vtable. A size no rule
/// set, as outside x86, is the real precision's: the 3M tile is its tile
/// (`params.zig`, where x86 takes the real `p` and `q` alike).
pub fn blocking(comptime T: type, vt: *const Vtable) Blocking {
    const real = Blocking.of(Float(T), vt);
    if (comptime portable.letter(T) == null) return real;
    const g = comptime portable.Gemm.of(T, true);
    const um: usize = @intCast(@max(1, portable.paramPtr(vt, g, .unroll_m).*));
    const un: usize = @intCast(@max(1, portable.paramPtr(vt, g, .unroll_n).*));
    const p = own(portable.paramPtr(vt, g, .p).*, real.p);
    return .{
        .p = @max(um, p / um * um),
        .q = own(portable.paramPtr(vt, g, .q).*, real.q),
        .r = @max(un, own(portable.paramPtr(vt, g, .r).*, real.r)),
        .um = um,
        .un = un,
    };
}

fn own(v: i32, fallback: usize) usize {
    return if (v > 0) @intCast(v) else fallback;
}

/// `gemm3m_level3.c`'s split of the remaining depth: all of it up to `max`,
/// half of it (rounded up) below twice `max`, else `max`.
pub fn depthSplit(rest: usize, max: usize) usize {
    if (rest >= 2 * max) return max;
    if (rest > max) return (rest + 1) / 2;
    return rest;
}

/// How the operands are packed (`driver.Pack`): `upper` and `lower` are the
/// symmetric operand's triangle, Hermitian with `herm`; `conj`: which `op`
/// conjugates (`gemm3m`'s `R` and `C`).
pub const Desc = struct {
    rows: Pack,
    cols: Pack,
    herm: bool = false,
    conj: driver.Conj = .n,
};

/// One pass: the part of `op(A)` and of `alpha * op(B)` it packs, and the
/// kernel's `alpha`, a sign per part of `C`.
const Pass = struct { a: Part, b: Part, re: i2, im: i2 };

/// The passes for each conjugation, `gemm3m_level3.c`'s `ALPHA5` to
/// `ALPHA18` and its choice of `OCOPYR` or `OCOPYI`; a conjugated `op(B)`
/// packs `conj(alpha) * B` (`-alpha[1]`).
fn passes(conj_a: bool, conj_b: bool) [3]Pass {
    if (conj_a == conj_b) {
        // `NN` and `RR`: the imaginary part from all three
        const s: i2 = if (conj_a) -1 else 1;
        return .{
            .{ .a = .b, .b = .b, .re = 0, .im = s },
            .{ .a = .r, .b = .r, .re = 1, .im = -s },
            .{ .a = .i, .b = .i, .re = -1, .im = -s },
        };
    }
    // `NR` and `RN`: the real part from all three
    const s: i2 = if (conj_a) 1 else -1;
    return .{
        .{ .a = .b, .b = .b, .re = 1, .im = 0 },
        .{ .a = .r, .b = .i, .re = -1, .im = s },
        .{ .a = .i, .b = .r, .re = -1, .im = -s },
    };
}

const raw = util.raw;

fn long(v: usize) BlasLong {
    return @intCast(v);
}

/// The kernel or copy `name` (an `Op`) of `T` on `args`.
inline fn call(comptime T: type, vt: *const Vtable, comptime name: []const u8, args: anytype) void {
    _ = @call(.auto, portable.kernel(T, vt, @field(Op, name)), args);
}

/// A symmetric copy, `pre` naming the routine and side.
inline fn symCall(comptime T: type, vt: *const Vtable, comptime pre: []const u8, upper: bool, comptime part: []const u8, args: anytype) void {
    if (upper) call(T, vt, pre ++ "ucopy" ++ part, args) else call(T, vt, pre ++ "lcopy" ++ part, args);
}

/// Packs part `part` of rows `is..` of `op(A)` over depth `ls..`
/// (`ICOPY*_OPERATION`).
fn packRows(comptime T: type, vt: *const Vtable, d: Desc, part: Part, x: Matrix(T), min_l: usize, min_i: usize, ls: usize, is: usize, dst: [*]Float(T)) void {
    switch (part) {
        inline else => |p| {
            const pt = @tagName(p);
            switch (d.rows) {
                .n => call(T, vt, "gemm3m_itcopy" ++ pt, .{ long(min_l), long(min_i), raw(T, x.p + is + ls * x.ld), long(x.ld), dst }),
                .t => call(T, vt, "gemm3m_incopy" ++ pt, .{ long(min_l), long(min_i), raw(T, x.p + ls + is * x.ld), long(x.ld), dst }),
                .upper, .lower => {
                    const args = .{ long(min_l), long(min_i), raw(T, x.p), long(x.ld), long(is), long(ls), dst };
                    if (d.herm) symCall(T, vt, "hemm3m_i", d.rows == .upper, pt, args) else symCall(T, vt, "symm3m_i", d.rows == .upper, pt, args);
                },
            }
        },
    }
}

/// Packs part `part` of `alpha` times columns `jjs..` of `op(B)` over depth
/// `ls..` (`OCOPY*_OPERATION`).
fn packCols(comptime T: type, vt: *const Vtable, d: Desc, part: Part, y: Matrix(T), min_l: usize, min_jj: usize, ls: usize, jjs: usize, ar: Float(T), ai: Float(T), dst: [*]Float(T)) void {
    switch (part) {
        inline else => |p| {
            const pt = @tagName(p);
            switch (d.cols) {
                .n => call(T, vt, "gemm3m_oncopy" ++ pt, .{ long(min_l), long(min_jj), raw(T, y.p + ls + jjs * y.ld), long(y.ld), ar, ai, dst }),
                .t => call(T, vt, "gemm3m_otcopy" ++ pt, .{ long(min_l), long(min_jj), raw(T, y.p + jjs + ls * y.ld), long(y.ld), ar, ai, dst }),
                .upper, .lower => {
                    const args = .{ long(min_l), long(min_jj), raw(T, y.p), long(y.ld), long(jjs), long(ls), ar, ai, dst };
                    if (d.herm) symCall(T, vt, "hemm3m_o", d.cols == .upper, pt, args) else symCall(T, vt, "symm3m_o", d.cols == .upper, pt, args);
                },
            }
        },
    }
}

/// The kernel on rows `[ib, ib + mi)` and columns `[jb, jb + nj)` of `C`
/// with `pass`'s signs, `C` scaled by `beta` first unless null.
fn tile(comptime T: type, vt: *const Vtable, pass: Pass, ib: usize, mi: usize, jb: usize, nj: usize, kl: usize, beta: ?T, sa: [*]Float(T), sb: [*]Float(T), c: [*]T, ldc: usize) void {
    const F = Float(T);
    const at = c + ib + jb * ldc;
    if (beta) |s| driver.scale(T, vt, mi, nj, s, at, ldc);
    const re: F = @floatFromInt(pass.re);
    const im: F = @floatFromInt(pass.im);
    call(T, vt, "gemm3m_kernel", .{ long(mi), long(nj), long(kl), re, im, sa, sb, raw(T, at), long(ldc) });
}

/// `C := beta*C + alpha * op(A) * op(B)` on rows `[m0, m1)` and columns
/// `[n0, n1)` of `C`, depth `k`: `gemm3m_level3.c` over `range_m` and
/// `range_n`. A null `beta` is 1.
pub fn drive(comptime T: type, vt: *const Vtable, d: Desc, blk: Blocking, m0: usize, m1: usize, n0: usize, n1: usize, k: usize, alpha: T, x: Matrix(T), y: Matrix(T), beta: ?T, c: [*]T, ldc: usize, work: driver.Work(Float(T))) void {
    const ps = passes(d.conj.ofA(), d.conj.ofB());
    const ai = if (d.conj.ofB()) -alpha.im else alpha.im;
    var js = n0;
    while (js < n1) : (js += blk.r) {
        const min_j = @min(blk.r, n1 - js);
        var ls: usize = 0;
        var min_l: usize = 0;
        while (ls < k) : (ls += min_l) {
            min_l = depthSplit(k - ls, blk.q);
            for (ps, 0..) |pass, pi| {
                const s = if (ls == 0 and pi == 0) beta else null;
                // the first row block packs `op(B)` as it goes
                var min_i = Blocking.split(m1 - m0, blk.p, blk.um);
                const once = min_i == m1 - m0;
                packRows(T, vt, d, pass.a, x, min_l, min_i, ls, m0, work.sa);
                var jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = @min(js + min_j - jjs, 3 * blk.un);
                    const pb = if (once) work.sb else work.sb + min_l * (jjs - js);
                    packCols(T, vt, d, pass.b, y, min_l, min_jj, ls, jjs, alpha.re, ai, pb);
                    tile(T, vt, pass, m0, min_i, jjs, min_jj, min_l, s, work.sa, pb, c, ldc);
                    jjs += min_jj;
                }
                var is = m0 + min_i;
                while (is < m1) : (is += min_i) {
                    min_i = Blocking.split(m1 - is, blk.p, blk.um);
                    packRows(T, vt, d, pass.a, x, min_l, min_i, ls, is, work.sa);
                    jjs = js;
                    while (jjs < js + min_j) {
                        const min_jj = @min(js + min_j - jjs, 3 * blk.un);
                        tile(T, vt, pass, is, min_i, jjs, min_jj, min_l, s, work.sa, work.sb + min_l * (jjs - js), c, ldc);
                        jjs += min_jj;
                    }
                }
            }
        }
    }
}

/// `drive` over all of `C` with packing space from `allocator` (a small
/// problem's on the stack); more than one worker `want`ed, over parts of `C`.
fn run(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, d: Desc, m: usize, n: usize, k: usize, alpha: T, x: Matrix(T), y: Matrix(T), beta: ?T, c: [*]T, ldc: usize, want: u32) Error!void {
    const F = Float(T);
    const W = driver.Work(F);
    const blk = blocking(T, vt);
    if (want > 1) return parts(allocator, io, T, vt, d, blk, m, n, k, alpha, x, y, beta, c, ldc, want);
    const lens = W.lens(blk, m, n, k, false);
    const len = W.total(lens);
    var stack: [driver.stack_bytes]u8 align(W.alignment) = undefined;
    const sc = try util.Scratch(F, W.alignment, W.alignment).init(allocator, &stack, len);
    defer sc.deinit(allocator);
    drive(T, vt, d, blk, 0, m, 0, n, k, alpha, x, y, beta, c, ldc, W.at(sc.ptr, lens));
}

/// `run` on up to `want` workers: the side of `C` with more panels split in
/// parts of whole panels (`unroll_n` columns or `unroll_m` rows, as
/// `gemm_thread_n` and `gemm_thread_m` divide it), each an independent
/// `drive` on packing space of its own, so its tiles are the serial ones.
fn parts(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, d: Desc, blk: Blocking, m: usize, n: usize, k: usize, alpha: T, x: Matrix(T), y: Matrix(T), beta: ?T, c: [*]T, ldc: usize, want: u32) Error!void {
    const F = Float(T);
    const W = driver.Work(F);
    const by_cols = (n + blk.un - 1) / blk.un >= (m + blk.um - 1) / blk.um;
    const Parts = struct {
        vt: *const Vtable,
        d: Desc,
        blk: Blocking,
        by_cols: bool,
        m: usize,
        n: usize,
        k: usize,
        alpha: T,
        x: Matrix(T),
        y: Matrix(T),
        beta: ?T,
        c: [*]T,
        ldc: usize,

        /// Part `ch`'s rows and columns of `C`.
        fn range(p: *const @This(), ch: threading.Chunk) [4]usize {
            return if (p.by_cols) .{ 0, p.m, ch.start, ch.end } else .{ ch.start, ch.end, 0, p.n };
        }

        fn lens(p: *const @This(), ch: threading.Chunk) [3]usize {
            const r = p.range(ch);
            return W.lens(p.blk, r[1] - r[0], r[3] - r[2], p.k, false);
        }

        fn len(p: *const @This(), ch: threading.Chunk) usize {
            return W.total(p.lens(ch));
        }

        fn one(p: *const @This(), ch: threading.Chunk, buf: []F) void {
            const r = p.range(ch);
            drive(T, p.vt, p.d, p.blk, r[0], r[1], r[2], r[3], p.k, p.alpha, p.x, p.y, p.beta, p.c, p.ldc, W.at(buf.ptr, p.lens(ch)));
        }
    };
    const plan: threading.Plan = .{
        .chunks = .init(if (by_cols) n else m, want, .even, if (by_cols) blk.un else blk.um),
        .workers = want,
    };
    const ps: Parts = .{ .vt = vt, .d = d, .blk = blk, .by_cols = by_cols, .m = m, .n = n, .k = k, .alpha = alpha, .x = x, .y = y, .beta = beta, .c = c, .ldc = ldc };
    try threading.forEachPart(allocator, io, &plan, F, W.alignment, &ps, Parts.len, Parts.one);
}

// ---- the routines ----------------------------------------------------------------

/// `C := alpha*op(A)*op(B) + beta*C` by 3M with the kernels of `vt`; `ta`,
/// `tb` transpose, `conj` conjugates (`interface/gemm.c` with `GEMM3M`).
/// Checked by the caller; `m`, `n` positive. Up to `want` workers.
pub fn gemm3m(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, ta: bool, tb: bool, conj: driver.Conj, m: usize, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize, want: u32) Error!void {
    const beta1 = driver.unlessOne(T, beta);
    if (k == 0 or util.isZero(T, alpha)) {
        if (beta1) |s| driver.scale(T, vt, m, n, s, c, ldc);
        return;
    }
    const x: Matrix(T) = .{ .p = @constCast(a), .ld = lda };
    const y: Matrix(T) = .{ .p = @constCast(b), .ld = ldb };
    const d: Desc = .{ .rows = if (ta) .t else .n, .cols = if (tb) .t else .n, .conj = conj };
    try run(allocator, io, T, vt, d, m, n, k, alpha, x, y, beta1, c, ldc, want);
}

/// `C := alpha*S*B + beta*C` (`left`) or `alpha*B*S + beta*C` by 3M, `S`
/// symmetric (`herm`: Hermitian) from its `upper` or lower triangle in `a`
/// (`symm3m_k.c`, `hemm3m_k.c`). `m`, `n` positive. Up to `want` workers.
pub fn symm3m(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, comptime herm: bool, left: bool, upper: bool, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize, want: u32) Error!void {
    const beta1 = driver.unlessOne(T, beta);
    if (util.isZero(T, alpha)) {
        if (beta1) |s| driver.scale(T, vt, m, n, s, c, ldc);
        return;
    }
    const s: Matrix(T) = .{ .p = @constCast(a), .ld = lda };
    const g: Matrix(T) = .{ .p = @constCast(b), .ld = ldb };
    const t: Pack = if (upper) .upper else .lower;
    const d: Desc = if (left) .{ .rows = t, .cols = .n, .herm = herm } else .{ .rows = .n, .cols = t, .herm = herm };
    try run(allocator, io, T, vt, d, m, n, if (left) m else n, alpha, if (left) s else g, if (left) g else s, beta1, c, ldc, want);
}

// ---- tests -----------------------------------------------------------------

const gemm_kernels = @import("gemm/kernel.zig");
const check = @import("../_test/check.zig");
const bitwise = @import("../_test/bitwise.zig");
const guard = @import("../_test/guard.zig");
const check3 = @import("../_test/blas3/check.zig");

test "3m passes: every conjugation adds op(A) times alpha op(B) once" {
    // per pass, the coefficient of each of ArBr, ArBi, AiBr, AiBi in each
    // part of `C`, by the parts the copies take (`b` both); `B` stands for
    // `alpha * B`, conjugated with `op(B)` as the driver packs it
    for ([_][2]bool{ .{ false, false }, .{ false, true }, .{ true, false }, .{ true, true } }) |cj| {
        var re: [2][2]i32 = @splat(@splat(0));
        var im: [2][2]i32 = @splat(@splat(0));
        for (passes(cj[0], cj[1])) |p| for (0..2) |x| for (0..2) |u| {
            const ax: i32 = if (p.a == .b or @backingInt(p.a) == x + 1) 1 else 0;
            const bu: i32 = if (p.b == .b or @backingInt(p.b) == u + 1) 1 else 0;
            re[x][u] += p.re * ax * bu;
            im[x][u] += p.im * ax * bu;
        };
        // `op(A) = Ar + sa i Ai`, `op(alpha B) = Br + sb i Bi`
        const sa: i32 = if (cj[0]) -1 else 1;
        const sb: i32 = if (cj[1]) -1 else 1;
        const want_re = [2][2]i32{ .{ 1, 0 }, .{ 0, -sa * sb } };
        const want_im = [2][2]i32{ .{ 0, sb }, .{ sa, 0 } };
        try std.testing.expectEqual(want_re, re);
        try std.testing.expectEqual(want_im, im);
    }
}

/// The portable kernels with the 3M blocking shrunk, so every block and
/// panel edge comes at small sizes.
fn shrunk(comptime T: type) Vtable {
    var vt: Vtable = .{};
    portable.install(&vt);
    const blk = blocking(T, &vt);
    check3.shrinkBlocking(T, &vt, true, blk.um, blk.un);
    return vt;
}

const fill = check.fillMixed;

const sameBits = bitwise.sameBits;

/// Each routine, shapes around every part and panel edge, at 1 to 9 workers
/// on `io` against one, bit for bit.
fn shapes(comptime T: type, io: Io) !void {
    const al = std.testing.allocator;
    const vt = shrunk(T);
    const blk = blocking(T, &vt);
    const dims = [_]usize{ 1, blk.um + 1, 2 * blk.p + blk.um + 1, blk.r + 1, 9 * blk.un + 5 };
    const ks = [_]usize{ 1, blk.q + 1, 2 * blk.q + 3 };
    const ld = 9 * blk.un + 2 * blk.p + blk.r + 2 * blk.q + 8;
    const a = try al.alloc(T, ld * ld);
    defer al.free(a);
    const b = try al.alloc(T, ld * ld);
    defer al.free(b);
    const c0 = try al.alloc(T, ld * ld);
    defer al.free(c0);
    const want = try al.alloc(T, c0.len);
    defer al.free(want);
    const got = try al.alloc(T, c0.len);
    defer al.free(got);
    fill(T, a, 1);
    fill(T, b, 2);
    fill(T, c0, 3);
    var poison: guard.Poison = .{ .child = al };
    const pa = poison.allocator();
    const alpha: T = .init(-0.75, 0.25);
    const beta: T = .init(0.5, -0.125);
    var x: usize = 0;
    for (0..8) |bits| for (dims) |m| for (dims) |n| for (ks) |k| {
        x += 1;
        if (x % 3 != 0) continue;
        const o: [3]bool = .{ bits & 1 != 0, bits & 2 != 0, bits & 4 != 0 };
        const Call = struct {
            fn f(w: u32, v: *const Vtable, al2: Allocator, io2: Io, oo: [3]bool, mm: usize, nn: usize, kk: usize, al_: T, be: T, ap: [*]const T, bp: [*]const T, l: usize, cp: [*]T) !void {
                // `gemm3m` over the transposes and conjugations, `symm3m` and
                // `hemm3m` over the sides and triangles
                switch (kk % 3) {
                    0 => try gemm3m(al2, io2, T, v, oo[0], oo[1], .of(oo[2], oo[0]), mm, nn, kk, al_, ap, l, bp, l, be, cp, l, w),
                    1 => try symm3m(al2, io2, T, v, false, oo[0], oo[1], mm, nn, al_, ap, l, bp, l, be, cp, l, w),
                    else => try symm3m(al2, io2, T, v, true, oo[0], oo[1], mm, nn, al_, ap, l, bp, l, be, cp, l, w),
                }
            }
        };
        @memcpy(want, c0);
        try Call.f(1, &vt, pa, Io.failing, o, m, n, k, alpha, beta, a.ptr, b.ptr, ld, want.ptr);
        for (2..10) |w| {
            @memcpy(got, c0);
            try Call.f(@intCast(w), &vt, pa, io, o, m, n, k, alpha, beta, a.ptr, b.ptr, ld, got.ptr);
            if (!sameBits(T, want, got)) {
                std.debug.print("{s} w={d} opt={any} m={d} n={d} k={d}\n", .{ @typeName(T), w, o, m, n, k });
                return error.TestUnexpectedResult;
            }
        }
    };
}

test "threaded 3m matches its serial bits at 1 to 9 workers" {
    inline for (.{ types.Complex32, types.Complex64 }) |T| {
        try shapes(T, Io.failing);
        if (comptime test_io.hosted) {
            var threaded: Io.Threaded = .init(std.testing.allocator, .{ .async_limit = .limited(3), .concurrent_limit = .limited(3) });
            defer threaded.deinit();
            try shapes(T, threaded.io());
        }
    }
}

test {
    _ = gemm_kernels;
}
