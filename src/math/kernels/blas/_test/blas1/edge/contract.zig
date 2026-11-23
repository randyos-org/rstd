//! One contract on every core: the extremes' NaN and index, and `nrm2`'s class, as the portable kernels give them.

const std = @import("std");
const types = @import("../../../types.zig");
const Vtable = @import("../../../Vtable.zig");
const portable = @import("../../../portable.zig");
const check = @import("../../check.zig");
const ed = @import("../edge.zig");

const BlasLong = types.BlasLong;
const C = std.math.Complex;
const Float = check.FloatOf;
const isComplex = check.isComplex;
const vtable_types = ed.vtable_types;

// ---- one contract on every core ----------------------------------------------

/// Long enough for every `@Vector` kernel's blocks (`4 * W` keys, 256 for
/// `f32` on AVX-512) to hold several.
pub const contract_n = 600;

/// Room for a `contract` operand.
pub var contract_mem: [contract_n]C(f64) = undefined;

/// Where `contract` puts special values in `x`: NaN at `nan`, and again at
/// `nan2`, infinities at `inf`, by element; null for none. `part` puts them
/// in the imaginary part of a complex element; `big` makes every element
/// with no special value `floatMax / 2` in each part, so a complex key
/// overflows.
pub const Specials = struct { nan: ?usize = null, nan2: ?usize = null, inf: ?usize = null, part: bool = false, big: bool = false };

/// Index and value routines (`i?amax`, `?amax`, ...) and `nrm2` against the
/// portable kernels on random data of every length to `contract_n`, with
/// NaN and infinities at the first, a middle and the last element, inside and
/// past the vector blocks: the same index, NaN or infinity where they give
/// one, and values within their accuracy class.
pub fn contract(vt: *const Vtable) !void {
    const ref = check.portableVtable();
    ed.failed = null;
    inline for (vtable_types) |p| check.keep(&ed.failed, contractOf(p[0], p[1], vt, ref));
    if (ed.failed) |e| return e;
}

pub fn contractOf(comptime l: []const u8, comptime T: type, vt: *const Vtable, ref: *const Vtable) !void {
    const F = Float(T);
    const complex = comptime isComplex(T);
    const w = if (complex) 2 else 1;
    const x: [*]F = @ptrCast(&contract_mem);
    for ([_]usize{ 1, 2, 9, 63, 257, 300, contract_n }) |n| {
        const mid = n / 2;
        const cases = [_]Specials{
            .{},
            .{ .nan = 0 },
            .{ .nan = mid },
            .{ .nan = n - 1 },
            .{ .nan = mid, .part = true },
            .{ .nan = mid, .nan2 = n - 1, .inf = 0 },
            .{ .inf = mid },
            .{ .inf = n - 1, .nan = n - 1, .part = true },
            .{ .big = true },
            .{ .big = true, .inf = n - 1, .part = true },
        };
        for (cases) |sp| {
            check.fill(F, x[0 .. w * n], n + 3);
            for (0..n) |i| {
                const at = w * i + @intFromBool(sp.part and complex);
                if (sp.big) for (0..w) |c| {
                    x[w * i + c] = std.math.copysign(std.math.floatMax(F) / 2, x[w * i + c]);
                };
                if (sp.inf == i) x[at] = -std.math.inf(F);
                if (sp.nan == i or sp.nan2 == i) x[at] = std.math.nan(F);
            }
            contractCase(l, T, vt, ref, n, x) catch |e| {
                check.say("{s} n={d} {any}: ", .{ @typeName(T), n, sp });
                return e;
            };
        }
    }
}

pub fn contractCase(comptime l: []const u8, comptime T: type, vt: *const Vtable, ref: *const Vtable, n: usize, x: [*]Float(T)) !void {
    const F = Float(T);
    const complex = comptime isComplex(T);
    const ni: BlasLong = @intCast(n);
    const prec: portable.Prec = comptime .from(l);
    inline for (check.extremes) |e| {
        if (comptime complex and !e[2].abs) continue;
        const i_got = portable.slotPtr(vt, prec, e[1]).*.?(ni, x, 1);
        const i_want = portable.slotPtr(ref, prec, e[1]).*.?(ni, x, 1);
        if (i_got != i_want) {
            check.say("{s}: {d}, the portable kernel {d}: ", .{ portable.slotName(prec, e[1]), i_got, i_want });
            return error.TestUnexpectedResult;
        }
        const got = portable.slotPtr(vt, prec, e[0]).*.?(ni, x, 1);
        const want = portable.slotPtr(ref, prec, e[0]).*.?(ni, x, 1);
        if (!(got == want or (std.math.isNan(got) and std.math.isNan(want)))) {
            check.say("{s}: {e}, the portable kernel {e}: ", .{ portable.slotName(prec, e[0]), got, want });
            return error.TestUnexpectedResult;
        }
    }
    const got = portable.slotPtr(vt, .from(l), .nrm2_k).*.?(ni, x, 1);
    const want = portable.slotPtr(ref, .from(l), .nrm2_k).*.?(ni, x, 1);
    const w: f128 = if (complex) 2 else 1;
    const ulps: f128 = if (F == f32) 2 else @as(f128, @floatFromInt(n)) * w + 4;
    if (std.math.isNan(want) != std.math.isNan(got) or std.math.isInf(want) != std.math.isInf(got) or
        (std.math.isFinite(want) and @abs(@as(f128, got) - want) > ulps * std.math.floatEps(F) * @as(f128, want)))
    {
        check.say("{s}nrm2: {e}, the portable kernel {e}: ", .{ l, got, want });
        return error.TestUnexpectedResult;
    }
}
