//! `SCHK1` to `SCHK5` of `test/[sdcz]blat3.f`: each group's sweep over the `.dat` values.

const std = @import("std");
const ob = @import("../__root__.zig");
const b3 = @import("../blat3.zig");
const mk = @import("make.zig");
const cl = @import("calls.zig");

const isComplex = ob.isComplex;
const testing = std.testing;
const Params = b3.Params;
const Routine = b3.Routine;
const nmax = b3.nmax;
const lseres = mk.lseres;
const make = mk.make;
const mmch = mk.mmch;
const mmch3m = mk.mmch3m;
const Call = cl.Call;
const Sweep = cl.Sweep;
const Work = cl.Work;
const bump = cl.bump;
const val = cl.val;

// ---- SCHK1 to SCHK5 -------------------------------------------------------------

/// `SCHK1`: `gemm` (and `gemm3m`).
pub fn chk1(comptime T: type, comptime r: Routine, p: Params) !void {
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    for (p.n) |m| for (p.n) |n| {
        const ldc = bump(m) orelse continue;
        const lcc = ldc * n;
        const null_ = n == 0 or m == 0;
        for (p.n) |k| for ("NTC") |ta| {
            const tra = ta != 'N';
            const ma = if (tra) k else m;
            const na = if (tra) m else k;
            const lda = bump(ma) orelse continue;
            const laa = lda * na;
            make(T, .ge, ' ', ' ', ma, na, w.a, nmax, w.aa, lda, &gen);
            for ("NTC") |tb| {
                const trb = tb != 'N';
                const mb = if (trb) n else k;
                const nb = if (trb) k else n;
                const ldb = bump(mb) orelse continue;
                const lbb = ldb * nb;
                make(T, .ge, ' ', ' ', mb, nb, w.b, nmax, w.bb, ldb, &gen);
                for (p.alpha) |alpha| for (p.beta) |beta| {
                    make(T, .ge, ' ', ' ', m, n, w.c, nmax, w.cc, ldc, &gen);
                    @memcpy(w.as[0..laa], w.aa[0..laa]);
                    @memcpy(w.bs[0..lbb], w.bb[0..lbb]);
                    @memcpy(w.cs[0..lcc], w.cc[0..lcc]);
                    const c: Call = .{ .transa = ta, .transb = tb, .m = m, .n = n, .k = k, .lda = lda, .ldb = ldb, .ldc = ldc, .alpha = alpha, .beta = beta };
                    try s.call(T, r, c, w);
                    try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                    try s.same(T, r, c, "B", ob.lse(T, w.bs[0..lbb], w.bb[0..lbb]));
                    try s.same(T, r, c, "C", if (null_) ob.lse(T, w.cs[0..lcc], w.cc[0..lcc]) else lseres(T, .ge, ' ', m, n, w.cc, w.cs, ldc));
                    if (null_) continue;
                    const mm = if (r == .gemm3m) mmch3m else mmch;
                    try s.result(T, r, c, mm(T, ta, tb, m, n, k, val(T, alpha), w.a, nmax, w.b, nmax, val(T, beta), w.c, nmax, w.ct, w.g, w.cc, ldc));
                };
            }
        };
    };
    try s.finish(T, r, p);
}

/// `SCHK2`: `symm` (complex: `symm`, `hemm`).
pub fn chk2(comptime T: type, comptime r: Routine, p: Params) !void {
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    for (p.n) |m| for (p.n) |n| {
        const ldc = bump(m) orelse continue;
        const lcc = ldc * n;
        const null_ = n == 0 or m == 0;
        const ldb = bump(m) orelse continue;
        const lbb = ldb * n;
        make(T, .ge, ' ', ' ', m, n, w.b, nmax, w.bb, ldb, &gen);
        for ("LR") |sd| {
            const left = sd == 'L';
            const na = if (left) m else n;
            const lda = bump(na) orelse continue;
            const laa = lda * na;
            for ("UL") |ul| {
                make(T, if (r == .hemm) .he else .sy, ul, ' ', na, na, w.a, nmax, w.aa, lda, &gen);
                for (p.alpha) |alpha| for (p.beta) |beta| {
                    make(T, .ge, ' ', ' ', m, n, w.c, nmax, w.cc, ldc, &gen);
                    @memcpy(w.as[0..laa], w.aa[0..laa]);
                    @memcpy(w.bs[0..lbb], w.bb[0..lbb]);
                    @memcpy(w.cs[0..lcc], w.cc[0..lcc]);
                    const c: Call = .{ .side = sd, .uplo = ul, .m = m, .n = n, .lda = lda, .ldb = ldb, .ldc = ldc, .alpha = alpha, .beta = beta };
                    try s.call(T, r, c, w);
                    try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                    try s.same(T, r, c, "B", ob.lse(T, w.bs[0..lbb], w.bb[0..lbb]));
                    try s.same(T, r, c, "C", if (null_) ob.lse(T, w.cs[0..lcc], w.cc[0..lcc]) else lseres(T, .ge, ' ', m, n, w.cc, w.cs, ldc));
                    if (null_) continue;
                    try s.result(T, r, c, if (left)
                        mmch(T, 'N', 'N', m, n, m, val(T, alpha), w.a, nmax, w.b, nmax, val(T, beta), w.c, nmax, w.ct, w.g, w.cc, ldc)
                    else
                        mmch(T, 'N', 'N', m, n, n, val(T, alpha), w.b, nmax, w.a, nmax, val(T, beta), w.c, nmax, w.ct, w.g, w.cc, ldc));
                };
            }
        }
    };
    try s.finish(T, r, p);
}

/// `SCHK3`: `trmm`, `trsm`.
pub fn chk3(comptime T: type, comptime r: Routine, p: Params) !void {
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    @memset(w.c, ob.zero(T));
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    const one = ob.val(T, 1, 0);
    const z = ob.zero(T);
    for (p.n) |m| next_n: for (p.n) |n| {
        const ldb = bump(m) orelse continue;
        const lbb = ldb * n;
        const null_ = m == 0 or n == 0;
        for ("LR") |sd| {
            const left = sd == 'L';
            const na = if (left) m else n;
            const lda = bump(na) orelse continue :next_n;
            const laa = lda * na;
            for ("UL") |ul| for ("NTC") |ta| for ("UN") |dg| for (p.alpha) |alpha| {
                make(T, .tr, ul, dg, na, na, w.a, nmax, w.aa, lda, &gen);
                make(T, .ge, ' ', ' ', m, n, w.b, nmax, w.bb, ldb, &gen);
                @memcpy(w.as[0..laa], w.aa[0..laa]);
                @memcpy(w.bs[0..lbb], w.bb[0..lbb]);
                const c: Call = .{ .side = sd, .uplo = ul, .transa = ta, .diag = dg, .m = m, .n = n, .lda = lda, .ldb = ldb, .alpha = alpha };
                try s.call(T, r, c, w);
                try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                try s.same(T, r, c, "B", if (null_) ob.lse(T, w.bs[0..lbb], w.bb[0..lbb]) else lseres(T, .ge, ' ', m, n, w.bb, w.bs, ldb));
                if (null_) continue;
                const av = val(T, alpha);
                if (r == .trmm) {
                    try s.result(T, r, c, if (left)
                        mmch(T, ta, 'N', m, n, m, av, w.a, nmax, w.b, nmax, z, w.c, nmax, w.ct, w.g, w.bb, ldb)
                    else
                        mmch(T, 'N', ta, m, n, n, av, w.b, nmax, w.a, nmax, z, w.c, nmax, w.ct, w.g, w.bb, ldb));
                } else {
                    // the computed solution into C, alpha*B where it was
                    for (0..n) |j| for (0..m) |i| {
                        w.c[i + j * nmax] = w.bb[i + j * ldb];
                        w.bb[i + j * ldb] = ob.mul(av, w.b[i + j * nmax]);
                    };
                    try s.result(T, r, c, if (left)
                        mmch(T, ta, 'N', m, n, m, one, w.a, nmax, w.c, nmax, z, w.b, nmax, w.ct, w.g, w.bb, ldb)
                    else
                        mmch(T, 'N', ta, m, n, n, one, w.c, nmax, w.a, nmax, z, w.b, nmax, w.ct, w.g, w.bb, ldb));
                }
            };
        }
    };
    try s.finish(T, r, p);
}

/// `SCHK4`: `syrk` (complex: `syrk` and `herk`, with real alpha and beta).
pub fn chk4(comptime T: type, comptime r: Routine, p: Params) !void {
    const herm = r == .herk;
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    const transes = if (comptime isComplex(T)) "NC" else "NTC";
    for (p.n) |n| {
        const ldc = bump(n) orelse continue;
        const lcc = ldc * n;
        for (p.n) |k| for (transes) |t0| {
            const tran = t0 != 'N';
            const tc: u8 = if (tran and comptime isComplex(T) and !herm) 'T' else t0;
            const ma = if (tran) k else n;
            const na = if (tran) n else k;
            const lda = bump(ma) orelse continue;
            const laa = lda * na;
            make(T, .ge, ' ', ' ', ma, na, w.a, nmax, w.aa, lda, &gen);
            for ("UL") |ul| {
                const up = ul == 'U';
                for (p.alpha) |alpha0| for (p.beta) |beta0| {
                    const alpha: [2]f64 = if (herm) .{ alpha0[0], 0 } else alpha0;
                    const beta: [2]f64 = if (herm) .{ beta0[0], 0 } else beta0;
                    var null_ = n == 0;
                    if (herm) null_ = null_ or ((k == 0 or alpha[0] == 0) and beta[0] == 1);
                    make(T, if (herm) .he else .sy, ul, ' ', n, n, w.c, nmax, w.cc, ldc, &gen);
                    @memcpy(w.as[0..laa], w.aa[0..laa]);
                    @memcpy(w.cs[0..lcc], w.cc[0..lcc]);
                    const c: Call = .{ .uplo = ul, .transa = tc, .n = n, .k = k, .lda = lda, .ldc = ldc, .alpha = alpha, .beta = beta };
                    try s.call(T, r, c, w);
                    try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                    try s.same(T, r, c, "C", if (null_) ob.lse(T, w.cs[0..lcc], w.cc[0..lcc]) else lseres(T, .sy, ul, n, n, w.cc, w.cs, ldc));
                    if (null_) continue;
                    const transt: u8 = if (herm) 'C' else 'T';
                    var jc: usize = 0;
                    for (0..n) |j| {
                        const jj = if (up) 0 else j;
                        const lj = if (up) j + 1 else n - j;
                        try s.result(T, r, c, if (tran)
                            mmch(T, transt, 'N', lj, 1, k, val(T, alpha), w.a[jj * nmax ..], nmax, w.a[j * nmax ..], nmax, val(T, beta), w.c[jj + j * nmax ..], nmax, w.ct, w.g, w.cc[jc..], ldc)
                        else
                            mmch(T, 'N', transt, lj, 1, k, val(T, alpha), w.a[jj..], nmax, w.a[j..], nmax, val(T, beta), w.c[jj + j * nmax ..], nmax, w.ct, w.g, w.cc[jc..], ldc));
                        jc += if (up) ldc else ldc + 1;
                    }
                };
            }
        };
    }
    try s.finish(T, r, p);
}

/// `SCHK5`: `syr2k` (complex: `syr2k` and `her2k`, with real beta). `ab`
/// holds `A` and `B` side by side (leading dimension `nmax`), or one over
/// the other (leading dimension `2*nmax`) when transposed.
pub fn chk5(comptime T: type, comptime r: Routine, p: Params) !void {
    const herm = r == .her2k;
    const cx = comptime isComplex(T);
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    const one = ob.val(T, 1, 0);
    const transes = if (cx) "NC" else "NTC";
    for (p.n) |n| {
        const ldc = bump(n) orelse continue;
        const lcc = ldc * n;
        for (p.n) |k| for (transes) |t0| {
            const tran = t0 != 'N';
            const tc: u8 = if (tran and cx and !herm) 'T' else t0;
            const ma = if (tran) k else n;
            const na = if (tran) n else k;
            const lda = bump(ma) orelse continue;
            const laa = lda * na;
            const ld_ab: usize = if (tran) 2 * nmax else nmax;
            make(T, .ge, ' ', ' ', ma, na, w.ab, ld_ab, w.aa, lda, &gen);
            const ldb = lda;
            const lbb = laa;
            make(T, .ge, ' ', ' ', ma, na, w.ab[if (tran) k else k * nmax..], ld_ab, w.bb, ldb, &gen);
            for ("UL") |ul| {
                const up = ul == 'U';
                for (p.alpha) |alpha| for (p.beta) |beta0| {
                    const beta: [2]f64 = if (herm) .{ beta0[0], 0 } else beta0;
                    const av = val(T, alpha);
                    var null_ = n == 0;
                    if (herm) null_ = null_ or ((k == 0 or ob.isZero(av)) and beta[0] == 1);
                    make(T, if (herm) .he else .sy, ul, ' ', n, n, w.c, nmax, w.cc, ldc, &gen);
                    @memcpy(w.as[0..laa], w.aa[0..laa]);
                    @memcpy(w.bs[0..lbb], w.bb[0..lbb]);
                    @memcpy(w.cs[0..lcc], w.cc[0..lcc]);
                    const c: Call = .{ .uplo = ul, .transa = tc, .n = n, .k = k, .lda = lda, .ldb = ldb, .ldc = ldc, .alpha = alpha, .beta = beta };
                    try s.call(T, r, c, w);
                    try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                    try s.same(T, r, c, "B", ob.lse(T, w.bs[0..lbb], w.bb[0..lbb]));
                    try s.same(T, r, c, "C", if (null_) ob.lse(T, w.cs[0..lcc], w.cc[0..lcc]) else lseres(T, .sy, ul, n, n, w.cc, w.cs, ldc));
                    if (null_) continue;
                    const transt: u8 = if (herm) 'C' else 'T';
                    // real `SCHK5` keeps alpha in the check; complex folds it into W
                    const ma_ = if (cx) one else av;
                    var jjab: usize = 0;
                    var jc: usize = 0;
                    for (0..n) |j| {
                        const jj = if (up) 0 else j;
                        const lj = if (up) j + 1 else n - j;
                        for (0..k) |i| {
                            const bv = if (tran) w.ab[j * 2 * nmax + k + i] else w.ab[(k + i) * nmax + j];
                            const avv = if (tran) w.ab[j * 2 * nmax + i] else w.ab[i * nmax + j];
                            if (!cx) {
                                w.w[i] = bv;
                                w.w[k + i] = avv;
                            } else if (tran) {
                                w.w[i] = ob.mul(av, bv);
                                w.w[k + i] = ob.mul(if (herm) ob.conj(av) else av, avv);
                            } else if (herm) {
                                w.w[i] = ob.mul(av, ob.conj(bv));
                                w.w[k + i] = ob.conj(ob.mul(av, avv));
                            } else {
                                w.w[i] = ob.mul(av, bv);
                                w.w[k + i] = ob.mul(av, avv);
                            }
                        }
                        try s.result(T, r, c, if (tran)
                            mmch(T, transt, 'N', lj, 1, 2 * k, ma_, w.ab[jjab..], 2 * nmax, w.w, 2 * nmax, val(T, beta), w.c[jj + j * nmax ..], nmax, w.ct, w.g, w.cc[jc..], ldc)
                        else
                            mmch(T, 'N', 'N', lj, 1, 2 * k, ma_, w.ab[jj..], nmax, w.w, 2 * nmax, val(T, beta), w.c[jj + j * nmax ..], nmax, w.ct, w.g, w.cc[jc..], ldc));
                        if (up) {
                            jc += ldc;
                        } else {
                            jc += ldc + 1;
                            if (tran) jjab += 2 * nmax;
                        }
                    }
                };
            }
        };
    }
    try s.finish(T, r, p);
}
