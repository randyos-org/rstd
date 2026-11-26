//! `zig run tools/probe/accel_quirks.zig -lc -framework Accelerate`, on a Mac: each Accelerate quirk
//! design.md lists ("Apple Accelerate", Quirks), called directly through each of Accelerate's three
//! symbol sets (`$NEWLAPACK$ILP64`, `$NEWLAPACK`, the classic names; `apple_accel.zig` binds the first
//! it finds), and whether this machine still shows it (issue #355). A zero stride on a written vector
//! makes Accelerate exit, so that case runs in a child.

const std = @import("std");
const sysinfo = @import("builtin");

const C = extern struct { re: f64, im: f64 };
const C32 = extern struct { re: f32, im: f32 };
const inf = std.math.inf(f64);
const nan = std.math.nan(f64);

// CBLAS's enums
const col = 102;
const no_trans = 111;
const upper = 121;
const lower = 122;
const non_unit = 131;
const left = 141;

/// One of Accelerate's symbol sets: its integer type and name suffix.
const Set = struct { I: type, suffix: []const u8, name: []const u8 };
const sets = .{
    Set{ .I = i64, .suffix = "$NEWLAPACK$ILP64", .name = "ILP64" },
    Set{ .I = c_int, .suffix = "$NEWLAPACK", .name = "NEWLAPACK" },
    Set{ .I = c_int, .suffix = "", .name = "classic" },
};

/// `name` in set `s`, or null where Accelerate lacks it.
fn sym(comptime s: Set, comptime F: type, comptime name: []const u8) ?*const F {
    return @extern(?*const F, .{ .name = name ++ s.suffix, .linkage = .weak });
}

/// LAPACK's names: `zpotrf_` classically, `zpotrf$NEWLAPACK` in the new sets.
fn lapack(comptime s: Set, comptime F: type, comptime base: []const u8) ?*const F {
    const name = if (s.suffix.len == 0) base ++ "_" else base ++ s.suffix;
    return @extern(?*const F, .{ .name = name, .linkage = .weak });
}

fn finite(x: f64) bool {
    return std.math.isFinite(x);
}

/// What a case found: the quirk, its absence, or no such symbol.
const Seen = enum { yes, no, absent };

fn of(b: bool) Seen {
    return if (b) .yes else .no;
}

fn Cases(comptime s: Set) type {
    const I = s.I;
    // past every kernel's blocks: vector paths differ from scalar ones
    const big = 67;
    const spot = 33;
    const big_dot = 100_000;
    return struct {
        /// The last case's detail, printed under its row.
        var note: []const u8 = "";
        var note_buf: [256]u8 = undefined;
        var first: []const u8 = "";
        var first_buf: [160]u8 = undefined;

        const list = [_]struct { []const u8, *const fn () Seen }{
            .{ "complex scal n=67, alpha (1.5, 0) on (inf, 1): imaginary part stays finite", zscal },
            .{ "complex axpy n=67, alpha (1.5, 0) on x (inf, 1): imaginary part stays finite", zaxpy },
            .{ "axpby, zero alpha: reads x (NaN reaches y)", daxpby },
            .{ "complex axpby n=67, alpha (1.5, 0) on x (inf, 1): imaginary part stays finite", zaxpby },
            .{ "complex gemv, beta (2, 0) on y (inf, 1): imaginary part stays finite", zgemv },
            .{ "trmv order 64, inf in x: rows below it get NaN", dtrmv },
            .{ "ger, inf in x, y(j) zero: column j gets NaN", dger },
            .{ "syr, inf in x, x(j) zero: column j gets NaN", dsyr },
            .{ "syr2, inf in x, x(j) and y(j) zero: column j gets NaN", dsyr2 },
            .{ "herk order 64, (inf, 1) in C, beta 0.5: imaginary part NaN", zherkBeta },
            .{ "herk: the diagonal keeps an imaginary part", zherkDiag },
            .{ "dsdot n=100000: off by more than f64 rounding", dsdot },
            .{ "gemm 67 x 64 x 64, inf in A's row 0: another row gets NaN", dgemm },
            .{ "zpotrf: a well conditioned matrix reported not positive definite", zpotrf },
            .{ "cpotrf: a well conditioned matrix reported not positive definite", cpotrf },
            .{ "trsm, subnormal diagonal: multiplies by its reciprocal", dtrsm },
            .{ "ssymv, inf in x, A nonzero: some y(i) NaN (issue #410)", ssymv },
            .{ "dsymv, inf in x, A nonzero: some y(i) NaN (issue #410)", dsymv },
        };

        fn zscal() Seen {
            const f = sym(s, fn (I, *const C, [*]C, I) callconv(.c) void, "cblas_zscal") orelse return .absent;
            var x: [big]C = @splat(.{ .re = 1, .im = 1 });
            x[spot] = .{ .re = inf, .im = 1 };
            f(big, &.{ .re = 1.5, .im = 0 }, &x, 1);
            return of(finite(x[spot].im));
        }

        fn zaxpy() Seen {
            const f = sym(s, fn (I, *const C, [*]const C, I, [*]C, I) callconv(.c) void, "cblas_zaxpy") orelse return .absent;
            var x: [big]C = @splat(.{ .re = 1, .im = 1 });
            x[spot] = .{ .re = inf, .im = 1 };
            var y: [big]C = @splat(.{ .re = 1, .im = 1 });
            f(big, &.{ .re = 1.5, .im = 0 }, &x, 1, &y, 1);
            return of(finite(y[spot].im));
        }

        fn daxpby() Seen {
            const f = sym(s, fn (I, f64, [*]const f64, I, f64, [*]f64, I) callconv(.c) void, "catlas_daxpby") orelse return .absent;
            const x = [_]f64{ nan, nan, nan, nan };
            var y = [_]f64{ 1, 2, 3, 4 };
            f(4, 0, &x, 1, 2, &y, 1);
            return of(std.math.isNan(y[0]));
        }

        fn zaxpby() Seen {
            const f = sym(s, fn (I, *const C, [*]const C, I, *const C, [*]C, I) callconv(.c) void, "catlas_zaxpby") orelse return .absent;
            var x: [big]C = @splat(.{ .re = 1, .im = 1 });
            x[spot] = .{ .re = inf, .im = 1 };
            var y: [big]C = @splat(.{ .re = 1, .im = 1 });
            f(big, &.{ .re = 1.5, .im = 0 }, &x, 1, &.{ .re = 1, .im = 1 }, &y, 1);
            return of(finite(y[spot].im));
        }

        fn zgemv() Seen {
            const f = sym(s, fn (c_int, c_int, I, I, *const C, [*]const C, I, [*]const C, I, *const C, [*]C, I) callconv(.c) void, "cblas_zgemv") orelse return .absent;
            const a = [_]C{.{ .re = 1, .im = 0 }};
            const x = [_]C{.{ .re = 1, .im = 0 }};
            var y = [_]C{.{ .re = inf, .im = 1 }};
            f(col, no_trans, 1, 1, &.{ .re = 1, .im = 0 }, &a, 1, &x, 1, &.{ .re = 2, .im = 0 }, &y, 1);
            return of(finite(y[0].im));
        }

        const n64 = 64;

        fn dtrmv() Seen {
            const f = sym(s, fn (c_int, c_int, c_int, c_int, I, [*]const f64, I, [*]f64, I) callconv(.c) void, "cblas_dtrmv") orelse return .absent;
            var a: [n64 * n64]f64 = @splat(1);
            var x: [n64]f64 = @splat(1);
            x[10] = inf;
            // upper, no transpose: row i reads x(j) for j >= i, so rows past 10 never meet it
            f(col, upper, no_trans, non_unit, n64, &a, n64, &x, 1);
            for (x[11..]) |v| if (!finite(v)) return .yes;
            return .no;
        }

        fn dger() Seen {
            const f = sym(s, fn (c_int, I, I, f64, [*]const f64, I, [*]const f64, I, [*]f64, I) callconv(.c) void, "cblas_dger") orelse return .absent;
            var a: [16]f64 = @splat(1);
            const x = [_]f64{ 1, inf, 1, 1 };
            const y = [_]f64{ 1, 0, 1, 1 };
            f(col, 4, 4, 1, &x, 1, &y, 1, &a, 4);
            for (a[4..8]) |v| if (!finite(v)) return .yes;
            return .no;
        }

        fn dsyr() Seen {
            const f = sym(s, fn (c_int, c_int, I, f64, [*]const f64, I, [*]f64, I) callconv(.c) void, "cblas_dsyr") orelse return .absent;
            var a: [16]f64 = @splat(1);
            const x = [_]f64{ 0, inf, 1, 1 };
            f(col, lower, 4, 1, &x, 1, &a, 4);
            for (a[0..4]) |v| if (!finite(v)) return .yes;
            return .no;
        }

        fn dsyr2() Seen {
            const f = sym(s, fn (c_int, c_int, I, f64, [*]const f64, I, [*]const f64, I, [*]f64, I) callconv(.c) void, "cblas_dsyr2") orelse return .absent;
            var a: [16]f64 = @splat(1);
            const x = [_]f64{ 0, inf, 1, 1 };
            const y = [_]f64{ 0, 1, 1, 1 };
            f(col, lower, 4, 1, &x, 1, &y, 1, &a, 4);
            for (a[0..4]) |v| if (!finite(v)) return .yes;
            return .no;
        }

        const Herk = fn (c_int, c_int, c_int, I, I, f64, [*]const C, I, f64, [*]C, I) callconv(.c) void;

        fn zherkBeta() Seen {
            const f = sym(s, Herk, "cblas_zherk") orelse return .absent;
            var a: [n64 * 4]C = @splat(.{ .re = 0.5, .im = 0.25 });
            var c: [n64 * n64]C = @splat(.{ .re = 1, .im = 0 });
            c[1] = .{ .re = inf, .im = 1 };
            f(col, lower, no_trans, n64, 4, 1, &a, n64, 0.5, &c, n64);
            return of(std.math.isNan(c[1].im));
        }

        fn zherkDiag() Seen {
            const f = sym(s, Herk, "cblas_zherk") orelse return .absent;
            var a: [3 * 4]C = @splat(.{ .re = 0.5, .im = 0.25 });
            var c: [9]C = @splat(.{ .re = 1, .im = 0 });
            for (0..3) |i| c[i * 4].im = 3;
            f(col, lower, no_trans, 3, 4, 1.5, &a, 3, 0.5, &c, 3);
            for (0..3) |i| if (c[i * 4].im != 0) return .yes;
            return .no;
        }

        fn dsdot() Seen {
            const f = sym(s, fn (I, [*]const f32, I, [*]const f32, I) callconv(.c) f64, "cblas_dsdot") orelse return .absent;
            var x: [big_dot]f32 = undefined;
            var y: [big_dot]f32 = undefined;
            var prng: std.Random.DefaultPrng = .init(7);
            const r = prng.random();
            for (&x, &y) |*a, *b| {
                a.* = r.float(f32) * 2 - 1;
                b.* = r.float(f32) * 2 - 1;
            }
            var want: f64 = 0;
            for (x, y) |a, b| want += @as(f64, a) * @as(f64, b);
            const got = f(big_dot, &x, 1, &y, 1);
            return of(@abs(got - want) > 1e-12 * @max(1, @abs(want)));
        }

        fn dgemm() Seen {
            const f = sym(s, fn (c_int, c_int, c_int, I, I, I, f64, [*]const f64, I, [*]const f64, I, f64, [*]f64, I) callconv(.c) void, "cblas_dgemm") orelse return .absent;
            const m = 67;
            const k = 64;
            const n = 64;
            var a: [m * k]f64 = @splat(1);
            a[0] = inf;
            const b: [k * n]f64 = @splat(1);
            var c: [m * n]f64 = @splat(0);
            f(col, no_trans, no_trans, m, n, k, 1, &a, m, &b, k, 0, &c, m);
            for (0..n) |j| for (1..m) |i| if (!finite(c[j * m + i])) return .yes;
            return .no;
        }

        fn zpotrf() Seen {
            const f = lapack(s, fn ([*]const u8, *const I, [*]C, *const I, *I) callconv(.c) void, "zpotrf") orelse return .absent;
            // 4 I plus ones: Hermitian, well conditioned
            var a: [16]C = @splat(.{ .re = 1, .im = 0 });
            for (0..4) |i| a[i * 5].re = 5;
            const n: I = 4;
            var info: I = 0;
            f("L", &n, &a, &n, &info);
            return of(info != 0);
        }

        fn cpotrf() Seen {
            const f = lapack(s, fn ([*]const u8, *const I, [*]C32, *const I, *I) callconv(.c) void, "cpotrf") orelse return .absent;
            var a: [16]C32 = @splat(.{ .re = 1, .im = 0 });
            for (0..4) |i| a[i * 5].re = 5;
            const n: I = 4;
            var info: I = 0;
            f("L", &n, &a, &n, &info);
            return of(info != 0);
        }

        fn dtrsm() Seen {
            const f = sym(s, fn (c_int, c_int, c_int, c_int, c_int, I, I, f64, [*]const f64, I, [*]f64, I) callconv(.c) void, "cblas_dtrsm") orelse return .absent;
            // A = diag(d), d subnormal; B = d, so A \ B is 1 by division, inf by 1/d
            const d = 1e-310;
            const a = [_]f64{d};
            var b = [_]f64{d};
            f(col, left, lower, no_trans, non_unit, 1, 1, 1, &a, 1, &b, 1);
            return of(b[0] != 1);
        }

        fn ssymv() Seen {
            return symv(f32, "cblas_ssymv");
        }

        fn dsymv() Seen {
            return symv(f64, "cblas_dsymv");
        }

        /// One infinity in `x` and no zero in `A`, so each `y(i)` is infinite;
        /// over orders, triangles, strides, the infinity's place and sign, and
        /// each buffer's offset (its alignment). `note` gets the failing count
        /// and the first failure.
        fn symv(comptime F: type, comptime name: []const u8) Seen {
            const f = sym(s, fn (c_int, c_int, I, F, [*]const F, I, [*]const F, I, F, [*]F, I) callconv(.c) void, name) orelse return .absent;
            const nmax = 130;
            var a: [nmax * nmax + 4]F = undefined;
            var x: [2 * nmax + 4]F = undefined;
            var y: [nmax + 4]F = undefined;
            var prng: std.Random.DefaultPrng = .init(11);
            const r = prng.random();
            for (&a) |*v| v.* = (r.float(F) + 0.5) * @as(F, if (r.boolean()) 1 else -1);
            for (&x) |*v| v.* = r.float(F) * 2 - 1;
            var bad: usize = 0;
            var runs: usize = 0;
            for ([_]usize{ 33, 130 }) |n| for ([_]c_int{ lower, upper }) |ul| for ([_]isize{ 1, -2 }) |incx| for ([_]usize{ 0, n / 2, n - 1 }) |p| for ([_]F{ std.math.inf(F), -std.math.inf(F) }) |v| for (0..4) |oa| for (0..4) |ox| for (0..4) |oy| {
                const step: usize = @abs(incx);
                const at = ox + step * (if (incx > 0) p else n - 1 - p);
                const keep = x[at];
                x[at] = v;
                @memset(&y, 0);
                f(col, ul, @intCast(n), 1, a[oa..].ptr, @intCast(n), x[ox..].ptr, @intCast(incx), 0, y[oy..].ptr, 1);
                x[at] = keep;
                runs += 1;
                const nan_at = for (y[oy..][0..n], 0..) |e, i| {
                    if (std.math.isNan(e)) break i;
                } else continue;
                if (bad == 0) first = std.mem.print(&first_buf, "n={d} {s} incx={d} inf at {d}, NaN at y({d}), offsets a {d} x {d} y {d}", .{ n, if (ul == lower) "lower" else "upper", incx, p, nan_at, oa, ox, oy }) catch "";
                bad += 1;
            };
            if (bad > 0) note = std.mem.print(&note_buf, "{d} of {d}; first {s}", .{ bad, runs, first }) catch first;
            return of(bad > 0);
        }

        /// The child's case: a zero stride on the vector `dswap` writes.
        fn zeroStride() void {
            const f = sym(s, fn (I, [*]f64, I, [*]f64, I) callconv(.c) void, "cblas_dswap") orelse std.process.exit(3);
            var x = [_]f64{ 1, 2, 3 };
            var y = [_]f64{ 4, 5, 6 };
            f(3, &x, 1, &y, 0);
        }
    };
}

pub fn main(init: std.process.Init) !void {
    if (comptime sysinfo.target.os.tag != .macos) @compileError("accel_quirks runs on macOS");
    const io = init.io;
    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len > 2 and std.mem.eql(u8, args[1], "zero-stride")) {
        inline for (sets) |s| if (std.mem.eql(u8, args[2], s.name)) Cases(s).zeroStride();
        return;
    }
    std.debug.print("{t}, Accelerate: quirk seen? (ILP64 | NEWLAPACK | classic)\n", .{sysinfo.target.cpu.arch});
    const n = Cases(sets[0]).list.len;
    for (0..n) |i| {
        std.debug.print("  ", .{});
        inline for (sets) |s| std.debug.print("{t:<6} ", .{Cases(s).list[i][1]()});
        std.debug.print("{s}\n", .{Cases(sets[0]).list[i][0]});
        inline for (sets) |s| if (Cases(s).note.len > 0) {
            std.debug.print("    {s}: {s}\n", .{ s.name, Cases(s).note });
            Cases(s).note = "";
        };
    }
    std.debug.print("  ", .{});
    inline for (sets) |s| {
        var child = try std.process.spawn(io, .{ .argv = &.{ args[0], "zero-stride", s.name }, .stdin = .ignore, .stdout = .ignore, .stderr = .ignore });
        const seen: Seen = switch (try child.wait(io)) {
            .exited => |code| if (code == 3) .absent else if (code != 0) .yes else .no,
            else => .yes,
        };
        std.debug.print("{t:<6} ", .{seen});
    }
    std.debug.print("zero stride on a written vector (swap): Accelerate exits\n", .{});
}
