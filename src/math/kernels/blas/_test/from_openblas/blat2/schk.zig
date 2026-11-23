//! `SCHK1` to `SCHK6` of `test/[sdcz]blat2.f`: each group's sweep over the `.dat` values.

const std = @import("std");
const ob = @import("../__root__.zig");
const b2 = @import("../blat2.zig");
const mk = @import("make.zig");
const cl = @import("calls.zig");

const isComplex = ob.isComplex;
const testing = std.testing;
const Params = b2.Params;
const Routine = b2.Routine;
const nmax = b2.nmax;
const lseres = mk.lseres;
const make = mk.make;
const makeVec = mk.makeVec;
const makeVecZeroMid = mk.makeVecZeroMid;
const mvch = mk.mvch;
const Call = cl.Call;
const Sweep = cl.Sweep;
const Work = cl.Work;

// ---- SCHK1 to SCHK6 -----------------------------------------------------------

/// `SCHK1`: `gemv` and `gbmv`.
pub fn chk1(comptime T: type, comptime r: Routine, p: Params) !void {
    const band = r == .gbmv;
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    for (p.n) |n| {
        const nd = n / 2 + 1;
        next_m: for (0..2) |im| {
            const m = if (im == 0) n -| nd else @min(n + nd, nmax);
            const nk = if (band) p.k.len else 1;
            for (0..nk) |iku| {
                const ku: isize = if (band) @intCast(p.k[iku]) else @as(isize, @intCast(n)) - 1;
                const kl: isize = if (band) @max(ku - 1, 0) else @as(isize, @intCast(m)) - 1;
                const lda = ob.bump(if (band) @intCast(kl + ku + 1) else m, nmax) orelse continue;
                const laa = lda * n;
                const null_ = n == 0 or m == 0;
                make(T, if (band) .gb else .ge, ' ', ' ', m, n, w.a, nmax, w.aa, lda, kl, ku, &gen, 0);
                for ("NTC") |trans| {
                    const tran = trans != 'N';
                    const ml = if (tran) n else m;
                    const nl = if (tran) m else n;
                    for (p.inc) |incx| {
                        const ax: usize = @abs(incx);
                        const lx = ax * nl;
                        makeVecZeroMid(T, nl, w.x, w.xx, ax, &gen, 0.5);
                        for (p.inc) |incy| {
                            const ay: usize = @abs(incy);
                            const ly = ay * ml;
                            for (p.alpha) |alpha| for (p.beta) |beta| {
                                makeVec(T, ml, w.y, w.yy, ay, &gen, 0);
                                @memcpy(w.as[0..laa], w.aa[0..laa]);
                                @memcpy(w.xs[0..lx], w.xx[0..lx]);
                                @memcpy(w.ys[0..ly], w.yy[0..ly]);
                                const c: Call = .{ .trans = trans, .m = m, .n = n, .k = kl, .ku = ku, .lda = lda, .incx = incx, .incy = incy, .alpha = alpha, .beta = beta };
                                try s.call(T, r, c, w);
                                try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                                try s.same(T, r, c, "X", ob.lse(T, w.xs[0..lx], w.xx[0..lx]));
                                try s.same(T, r, c, "Y", if (null_) ob.lse(T, w.ys[0..ly], w.yy[0..ly]) else lseres(T, .ge, ' ', 1, ml, w.yy, w.ys, ay));
                                if (null_) continue :next_m;
                                try s.result(T, r, c, mvch(T, trans, m, n, ob.val(T, alpha[0], alpha[1]), w.a, nmax, w.x, incx, ob.val(T, beta[0], beta[1]), w.y, incy, w.yt, w.g, w.yy));
                            };
                        }
                    }
                }
            }
        }
    }
    try s.finish(T, r, p);
}

/// `SCHK2`: `symv`, `sbmv`, `spmv` (complex: `hemv`, `hbmv`, `hpmv`).
pub fn chk2(comptime T: type, comptime r: Routine, p: Params) !void {
    const band = r == .sbmv;
    const pack = r == .spmv;
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    next_n: for (p.n) |n| {
        const nk = if (band) p.k.len else 1;
        for (0..nk) |ik| {
            const k: isize = if (band) @intCast(p.k[ik]) else @as(isize, @intCast(n)) - 1;
            const lda = ob.bump(if (band) @intCast(k + 1) else n, nmax) orelse continue;
            const laa = if (pack) n * (n + 1) / 2 else lda * n;
            const null_ = n == 0;
            for ("UL") |uplo| {
                make(T, if (band) .sb else if (pack) .sp else .sy, uplo, ' ', n, n, w.a, nmax, w.aa, lda, k, k, &gen, 0);
                for (p.inc) |incx| {
                    const ax: usize = @abs(incx);
                    const lx = ax * n;
                    makeVecZeroMid(T, n, w.x, w.xx, ax, &gen, 0.5);
                    for (p.inc) |incy| {
                        const ay: usize = @abs(incy);
                        const ly = ay * n;
                        for (p.alpha) |alpha| for (p.beta) |beta| {
                            makeVec(T, n, w.y, w.yy, ay, &gen, 0);
                            @memcpy(w.as[0..laa], w.aa[0..laa]);
                            @memcpy(w.xs[0..lx], w.xx[0..lx]);
                            @memcpy(w.ys[0..ly], w.yy[0..ly]);
                            const c: Call = .{ .uplo = uplo, .n = n, .k = k, .lda = lda, .incx = incx, .incy = incy, .alpha = alpha, .beta = beta };
                            try s.call(T, r, c, w);
                            try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                            try s.same(T, r, c, "X", ob.lse(T, w.xs[0..lx], w.xx[0..lx]));
                            try s.same(T, r, c, "Y", if (null_) ob.lse(T, w.ys[0..ly], w.yy[0..ly]) else lseres(T, .ge, ' ', 1, n, w.yy, w.ys, ay));
                            if (null_) continue :next_n;
                            try s.result(T, r, c, mvch(T, 'N', n, n, ob.val(T, alpha[0], alpha[1]), w.a, nmax, w.x, incx, ob.val(T, beta[0], beta[1]), w.y, incy, w.yt, w.g, w.yy));
                        };
                    }
                }
            }
        }
    }
    try s.finish(T, r, p);
}

/// `SCHK3`: `trmv`, `tbmv`, `tpmv`, `trsv`, `tbsv`, `tpsv`.
pub fn chk3(comptime T: type, comptime r: Routine, p: Params) !void {
    const band = r == .tbmv or r == .tbsv;
    const pack = r == .tpmv or r == .tpsv;
    const solve = r == .trsv or r == .tbsv or r == .tpsv;
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    @memset(w.z, ob.zero(T));
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    const one = ob.val(T, 1, 0);
    next_n: for (p.n) |n| {
        const nk = if (band) p.k.len else 1;
        for (0..nk) |ik| {
            const k: isize = if (band) @intCast(p.k[ik]) else @as(isize, @intCast(n)) - 1;
            const lda = ob.bump(if (band) @intCast(k + 1) else n, nmax) orelse continue;
            const laa = if (pack) n * (n + 1) / 2 else lda * n;
            const null_ = n == 0;
            for ("UL") |uplo| for ("NTC") |trans| for ("UN") |diag| {
                make(T, if (band) .tb else if (pack) .tp else .tr, uplo, diag, n, n, w.a, nmax, w.aa, lda, k, k, &gen, 0);
                for (p.inc) |incx| {
                    const ax: usize = @abs(incx);
                    const lx = ax * n;
                    makeVecZeroMid(T, n, w.x, w.xx, ax, &gen, 0.5);
                    @memcpy(w.as[0..laa], w.aa[0..laa]);
                    @memcpy(w.xs[0..lx], w.xx[0..lx]);
                    const c: Call = .{ .uplo = uplo, .trans = trans, .diag = diag, .n = n, .k = k, .lda = lda, .incx = incx };
                    try s.call(T, r, c, w);
                    try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                    try s.same(T, r, c, "X", if (null_) ob.lse(T, w.xs[0..lx], w.xx[0..lx]) else lseres(T, .ge, ' ', 1, n, w.xx, w.xs, ax));
                    if (null_) continue :next_n;
                    if (!solve) {
                        try s.result(T, r, c, mvch(T, trans, n, n, one, w.a, nmax, w.x, incx, ob.zero(T), w.z, incx, w.yt, w.g, w.xx));
                    } else {
                        // the approximation to the original vector
                        for (0..n) |i| {
                            w.z[i] = w.xx[i * ax];
                            w.xx[i * ax] = w.x[i];
                        }
                        try s.result(T, r, c, mvch(T, trans, n, n, one, w.a, nmax, w.z, incx, ob.zero(T), w.x, incx, w.yt, w.g, w.xx));
                    }
                }
            };
        }
    }
    try s.finish(T, r, p);
}

/// `SCHK4`: `ger` (complex: `geru`, `gerc`).
pub fn chk4(comptime T: type, comptime r: Routine, p: Params) !void {
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    const one = ob.val(T, 1, 0);
    for (p.n) |n| {
        const nd = n / 2 + 1;
        next_m: for (0..2) |im| {
            const m = if (im == 0) n -| nd else @min(n + nd, nmax);
            const lda = ob.bump(m, nmax) orelse continue;
            const laa = lda * n;
            const null_ = n == 0 or m == 0;
            for (p.inc) |incx| {
                const ax: usize = @abs(incx);
                const lx = ax * m;
                makeVecZeroMid(T, m, w.x, w.xx, ax, &gen, 0.5);
                for (p.inc) |incy| {
                    const ay: usize = @abs(incy);
                    const ly = ay * n;
                    makeVecZeroMid(T, n, w.y, w.yy, ay, &gen, 0);
                    for (p.alpha) |alpha| {
                        make(T, .ge, ' ', ' ', m, n, w.a, nmax, w.aa, lda, @as(isize, @intCast(m)) - 1, @as(isize, @intCast(n)) - 1, &gen, 0);
                        @memcpy(w.as[0..laa], w.aa[0..laa]);
                        @memcpy(w.xs[0..lx], w.xx[0..lx]);
                        @memcpy(w.ys[0..ly], w.yy[0..ly]);
                        const c: Call = .{ .m = m, .n = n, .lda = lda, .incx = incx, .incy = incy, .alpha = alpha };
                        try s.call(T, r, c, w);
                        try s.same(T, r, c, "X", ob.lse(T, w.xs[0..lx], w.xx[0..lx]));
                        try s.same(T, r, c, "Y", ob.lse(T, w.ys[0..ly], w.yy[0..ly]));
                        try s.same(T, r, c, "A", if (null_) ob.lse(T, w.as[0..laa], w.aa[0..laa]) else lseres(T, .ge, ' ', m, n, w.aa, w.as, lda));
                        if (null_) continue :next_m;
                        const z = w.z[0..m];
                        for (z, 0..) |*e, i| e.* = w.x[if (incx > 0) i else m - 1 - i];
                        for (0..n) |j| {
                            var wv = [1]T{w.y[if (incy > 0) j else n - 1 - j]};
                            if (r == .gerc) wv[0] = ob.conj(wv[0]);
                            try s.result(T, r, c, mvch(T, 'N', m, 1, ob.val(T, alpha[0], alpha[1]), z, nmax, &wv, 1, one, w.a[j * nmax ..], 1, w.yt, w.g, w.aa[j * lda ..]));
                        }
                    }
                }
            }
        }
    }
    try s.finish(T, r, p);
}

/// `SCHK5`: `syr`, `spr` (complex: `her`, `hpr`, with the real part of
/// alpha).
pub fn chk5(comptime T: type, comptime r: Routine, p: Params) !void {
    const pack = r == .spr;
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    const one = ob.val(T, 1, 0);
    next_n: for (p.n) |n| {
        const lda = ob.bump(n, nmax) orelse continue;
        const laa = if (pack) n * (n + 1) / 2 else lda * n;
        for ("UL") |uplo| {
            const up = uplo == 'U';
            for (p.inc) |incx| {
                const ax: usize = @abs(incx);
                const lx = ax * n;
                makeVecZeroMid(T, n, w.x, w.xx, ax, &gen, 0.5);
                for (p.alpha) |alpha_c| {
                    const alpha: [2]f64 = .{ alpha_c[0], 0 };
                    const null_ = n == 0 or alpha[0] == 0;
                    make(T, if (pack) .sp else .sy, uplo, ' ', n, n, w.a, nmax, w.aa, lda, @as(isize, @intCast(n)) - 1, @as(isize, @intCast(n)) - 1, &gen, 0);
                    @memcpy(w.as[0..laa], w.aa[0..laa]);
                    @memcpy(w.xs[0..lx], w.xx[0..lx]);
                    const c: Call = .{ .uplo = uplo, .n = n, .lda = lda, .incx = incx, .alpha = alpha };
                    try s.call(T, r, c, w);
                    try s.same(T, r, c, "X", ob.lse(T, w.xs[0..lx], w.xx[0..lx]));
                    try s.same(T, r, c, "A", if (null_) ob.lse(T, w.as[0..laa], w.aa[0..laa]) else pack or lseres(T, .sy, uplo, n, n, w.aa, w.as, lda));
                    if (null_) {
                        if (n == 0) continue :next_n;
                        continue;
                    }
                    const z = w.z[0..n];
                    for (z, 0..) |*e, i| e.* = w.x[if (incx > 0) i else n - 1 - i];
                    var ja: usize = 0;
                    for (0..n) |j| {
                        const wv = [1]T{ob.conj(z[j])};
                        const jj = if (up) 0 else j;
                        const lj = if (up) j + 1 else n - j;
                        try s.result(T, r, c, mvch(T, 'N', lj, 1, ob.val(T, alpha[0], alpha[1]), z[jj..], lj, &wv, 1, one, w.a[jj + j * nmax ..], 1, w.yt, w.g, w.aa[ja..]));
                        ja += if (pack) lj else if (up) lda else lda + 1;
                    }
                }
            }
        }
    }
    try s.finish(T, r, p);
}

/// `SCHK6`: `syr2`, `spr2` (complex: `her2`, `hpr2`).
pub fn chk6(comptime T: type, comptime r: Routine, p: Params) !void {
    const pack = r == .spr2;
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    const one = ob.val(T, 1, 0);
    next_n: for (p.n) |n| {
        const lda = ob.bump(n, nmax) orelse continue;
        const laa = if (pack) n * (n + 1) / 2 else lda * n;
        for ("UL") |uplo| {
            const up = uplo == 'U';
            for (p.inc) |incx| {
                const ax: usize = @abs(incx);
                const lx = ax * n;
                makeVecZeroMid(T, n, w.x, w.xx, ax, &gen, 0.5);
                for (p.inc) |incy| {
                    const ay: usize = @abs(incy);
                    const ly = ay * n;
                    makeVecZeroMid(T, n, w.y, w.yy, ay, &gen, 0);
                    for (p.alpha) |alpha| {
                        const av = ob.val(T, alpha[0], alpha[1]);
                        const null_ = n == 0 or ob.isZero(av);
                        make(T, if (pack) .sp else .sy, uplo, ' ', n, n, w.a, nmax, w.aa, lda, @as(isize, @intCast(n)) - 1, @as(isize, @intCast(n)) - 1, &gen, 0);
                        @memcpy(w.as[0..laa], w.aa[0..laa]);
                        @memcpy(w.xs[0..lx], w.xx[0..lx]);
                        @memcpy(w.ys[0..ly], w.yy[0..ly]);
                        const c: Call = .{ .uplo = uplo, .n = n, .lda = lda, .incx = incx, .incy = incy, .alpha = alpha };
                        try s.call(T, r, c, w);
                        try s.same(T, r, c, "X", ob.lse(T, w.xs[0..lx], w.xx[0..lx]));
                        try s.same(T, r, c, "Y", ob.lse(T, w.ys[0..ly], w.yy[0..ly]));
                        try s.same(T, r, c, "A", if (null_) ob.lse(T, w.as[0..laa], w.aa[0..laa]) else pack or lseres(T, .sy, uplo, n, n, w.aa, w.as, lda));
                        if (null_) {
                            if (n == 0) continue :next_n;
                            continue;
                        }
                        // Z(:, 1) is x, Z(:, 2) is y, both forwards
                        const z1 = w.z[0..n];
                        const z2 = w.z[nmax..][0..n];
                        for (0..n) |i| {
                            z1[i] = w.x[if (incx > 0) i else n - 1 - i];
                            z2[i] = w.y[if (incy > 0) i else n - 1 - i];
                        }
                        var ja: usize = 0;
                        for (0..n) |j| {
                            const wv: [2]T = if (comptime isComplex(T))
                                .{ ob.mul(av, ob.conj(z2[j])), ob.mul(ob.conj(av), ob.conj(z1[j])) }
                            else
                                .{ z2[j], z1[j] };
                            const ma = if (comptime isComplex(T)) one else av;
                            const jj = if (up) 0 else j;
                            const lj = if (up) j + 1 else n - j;
                            try s.result(T, r, c, mvch(T, 'N', lj, 2, ma, w.z[jj..], nmax, &wv, 1, one, w.a[jj + j * nmax ..], 1, w.yt, w.g, w.aa[ja..]));
                            ja += if (pack) lj else if (up) lda else lda + 1;
                        }
                    }
                }
            }
        }
    }
    try s.finish(T, r, p);
}
