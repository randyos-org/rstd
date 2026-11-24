//! The edge sweep's operands: each routine's matrix kind, layout and shape, and the matrix that holds them.

const std = @import("std");
const check = @import("../../check.zig");
const bitwise = @import("../../bitwise.zig");
const pattern = @import("../../pattern.zig");
const reference = @import("../../reference.zig");

const Allocator = std.mem.Allocator;
const Float = check.FloatOf;
const isComplex = check.isComplex;
const Wide = reference.Wide;
const Pattern = pattern.Pattern;

pub const Kind = enum { general, symmetric, hermitian, hermitian_rev, triangular };
pub const Layout = enum { full, band, packed_ };

/// What a matrix operand is and how it is stored.
pub const Shape = struct {
    kind: Kind,
    layout: Layout = .full,
    m: usize,
    n: usize,
    /// Sub- and superdiagonals of a band; a triangle's other one is 0.
    kl: usize = 0,
    ku: usize = 0,
    lower: bool = false,
    unit: bool = false,
    /// A solve's: the diagonal dominant and left out of the pattern.
    solve: bool = false,
    /// Rows of padding below a full or band matrix.
    pad: usize = 1,

    pub fn ld(s: Shape) usize {
        return switch (s.layout) {
            .full => s.m + s.pad,
            .band => s.kl + s.ku + 1 + s.pad,
            .packed_ => 1,
        };
    }

    pub fn len(s: Shape) usize {
        return if (s.layout == .packed_) s.n * (s.n + 1) / 2 else s.ld() * s.n;
    }

    /// Where `(i, j)` is stored, if it is.
    pub fn stored(s: Shape, i: usize, j: usize) ?usize {
        if (s.kind == .general) return switch (s.layout) {
            .full => i + j * s.ld(),
            .band => if (i + s.ku >= j and i <= j + s.kl) s.ku + i - j + j * s.ld() else null,
            .packed_ => unreachable,
        };
        if (if (s.lower) i < j else i > j) return null;
        const d = if (s.lower) i - j else j - i;
        return switch (s.layout) {
            .full => i + j * s.ld(),
            .band => if (d > s.kl + s.ku) null else (if (s.lower) d else s.ku - d) + j * s.ld(),
            .packed_ => if (s.lower) j * (2 * s.n - j + 1) / 2 + d else j * (j + 1) / 2 + i,
        };
    }

    /// Whether the routine reads entry `(i, j)`: stored, and not a unit diagonal.
    pub fn read(s: Shape, i: usize, j: usize) ?usize {
        if (i == j and s.unit) return null;
        return s.stored(i, j);
    }

    /// The entries a pattern goes through, in `init`'s order.
    pub fn entries(s: Shape) usize {
        var k: usize = 0;
        for (0..s.n) |j| for (0..s.m) |i| {
            if (s.read(i, j) != null and !(s.solve and i == j)) k += 1;
        };
        return k;
    }
};

pub fn Matrix(comptime T: type) type {
    return struct {
        const Self = @This();
        s: Shape,
        buf: []T,
        was: []T,

        /// Random entries where the routine reads, shaped by `p` in order;
        /// NaN everywhere else, and in a Hermitian diagonal's imaginary part.
        pub fn init(al: Allocator, s: Shape, seed: u64, p: ?Pattern) !Self {
            const F = Float(T);
            const buf = try al.alloc(T, s.len() + check.tail);
            errdefer al.free(buf);
            @memset(std.mem.sliceAsBytes(buf), 0xff);
            const vals = try al.alloc(T, s.len());
            defer al.free(vals);
            check.fill(T, vals, seed);
            const sc: F = if (s.solve) 1 / @sqrt(@as(F, @floatFromInt(s.n + 1))) else 1;
            var k: usize = 0;
            for (0..s.n) |j| for (0..s.m) |i| {
                const idx = s.read(i, j) orelse continue;
                var v = vals[idx];
                if (s.solve and i == j) {
                    v = if (comptime isComplex(T)) .init(3 + @abs(v.re), v.im) else 3 + @abs(v);
                } else {
                    if (s.solve) v = if (comptime isComplex(T)) .init(v.re * sc, v.im * sc) else v * sc;
                    if (p) |q| v = pattern.shape(T, v, k, q);
                    k += 1;
                }
                if (comptime isComplex(T)) if (i == j and (s.kind == .hermitian or s.kind == .hermitian_rev)) {
                    buf[idx].re = v.re;
                    continue;
                };
                buf[idx] = v;
            };
            return .{ .s = s, .buf = buf, .was = try al.dupe(T, buf) };
        }

        pub fn deinit(self: *Self, al: Allocator) void {
            al.free(self.was);
            al.free(self.buf);
        }

        pub fn ld(self: *const Self) isize {
            return @intCast(self.s.ld());
        }

        /// NaN in every entry: the routine must not read them.
        pub fn poison(self: *Self) void {
            @memset(std.mem.sliceAsBytes(self.buf), 0xff);
            @memcpy(self.was, self.buf);
        }

        /// Entry `(i, j)` of the matrix the storage stands for, as it was;
        /// null for a zero outside the band or triangle, which the reference
        /// never multiplies.
        pub fn at(self: *const Self, i: usize, j: usize) ?Wide(T) {
            const s = self.s;
            const z: ?Wide(T) = null;
            switch (s.kind) {
                .general => return if (s.stored(i, j)) |k| reference.wide(T, self.was[k]) else z,
                .triangular => {
                    const k = s.stored(i, j) orelse return z;
                    if (i == j and s.unit) return reference.wide(T, if (comptime isComplex(T)) T.init(1, 0) else 1);
                    return reference.wide(T, self.was[k]);
                },
                .symmetric => {
                    const k = s.stored(i, j) orelse s.stored(j, i) orelse return z;
                    return reference.wide(T, self.was[k]);
                },
                .hermitian, .hermitian_rev => {
                    if (comptime !isComplex(T)) unreachable;
                    if (i == j) return .init(self.was[s.stored(i, i).?].re, 0);
                    const rev = s.kind == .hermitian_rev;
                    if (s.stored(i, j)) |k| return if (rev) reference.conjW(reference.wide(T, self.was[k])) else reference.wide(T, self.was[k]);
                    const k = s.stored(j, i) orelse return z;
                    return if (rev) reference.wide(T, self.was[k]) else reference.conjW(reference.wide(T, self.was[k]));
                },
            }
        }

        pub fn unchanged(self: *const Self) !void {
            for (self.buf, self.was, 0..) |a, b, k| if (!bitwise.same(T, a, b, true)) {
                check.say("matrix storage {d} changed: ", .{k});
                return error.TestUnexpectedResult;
            };
        }
    };
}
