//! The same bits at any address: each kernel on operands shifted through every misalignment.

const std = @import("std");
const types = @import("../../../types.zig");
const Vtable = @import("../../../Vtable.zig");
const portable = @import("../../../portable.zig");
const util = @import("../../../util.zig");
const check = @import("../../check.zig");
const bitwise = @import("../../bitwise.zig");
const ed = @import("../edge.zig");
const er = @import("runs.zig");

const BlasLong = types.BlasLong;
const C = std.math.Complex;
const Float = check.FloatOf;
const isComplex = check.isComplex;
const vtable_types = ed.vtable_types;
const Reduce = er.Reduce;

// ---- the same bits at any address ----------------------------------------------

/// Room for both `shifted` operands at the largest shift, 64-byte aligned so
/// each shift is the same misalignment on every run.
pub var shift_mem: [2 * (offset_n + 8)]C(f64) align(64) = undefined;
pub const offset_n = 300;

/// `kern` (a reduction of one or two vectors, by `two`) on the same data
/// shifted by 0 to 7 floats in each operand, any alignment a float may have:
/// the same bits each time, as `threading.md`'s reproducibility needs.
pub fn shifted(comptime T: type, routine: []const u8, comptime two: bool, kern: anytype) !void {
    const F = Float(T);
    const mem: [*]F = @ptrCast(&shift_mem);
    const w = @sizeOf(T) / @sizeOf(F);
    var vals: [2][offset_n]T = undefined;
    check.fill(T, &vals[0], 7);
    check.fill(T, &vals[1], 8);
    for ([_]usize{ 1, 3, 7, 16, 33, 67, 100, offset_n }) |n| {
        var first: ?[2]f64 = null;
        for (0..8) |ox| for (0..if (two) 8 else 1) |oy| {
            const x = mem + ox;
            const y = mem + w * (offset_n + 8) + oy;
            @memcpy(x[0 .. w * n], @as([*]const F, @ptrCast(&vals[0]))[0 .. w * n]);
            @memcpy(y[0 .. w * n], @as([*]const F, @ptrCast(&vals[1]))[0 .. w * n]);
            const got = if (two) kern(@intCast(n), x, 1, y, 1) else kern(@intCast(n), x, 1);
            const bits: [2]f64 = if (@typeInfo(@TypeOf(got)) == .@"struct") .{ got.re, got.im } else .{ got, 0 };
            if (first) |f| {
                if (!bitwise.same(f64, f[0], bits[0], true) or !bitwise.same(f64, f[1], bits[1], true)) {
                    check.say("{s} {s} n={d}: shifted by {d},{d} gave {any}, unshifted {any}: ", .{ routine, @typeName(T), n, ox, oy, bits, f });
                    return error.TestUnexpectedResult;
                }
            } else first = bits;
        };
    }
}

/// The element-wise kernels `shiftedWise` runs.
pub const Wise = enum { scal, axpy, axpby, rot };

/// `kern` (an element-wise kernel, by `op`, with `alpha` `a`) on `shifted`'s
/// data and shifts: the same bits in `x` and `y` each time.
pub fn shiftedWise(comptime T: type, routine: []const u8, comptime op: Wise, kern: anytype, a: T) !void {
    const F = Float(T);
    const complex = comptime isComplex(T);
    const mem: [*]F = @ptrCast(&shift_mem);
    const w = @sizeOf(T) / @sizeOf(F);
    var vals: [2][offset_n]T = undefined;
    check.fill(T, &vals[0], 7);
    check.fill(T, &vals[1], 8);
    const b: T = if (complex) .init(-0.5, 0.25) else -0.5;
    var first: [2][w * offset_n]F = undefined;
    for ([_]usize{ 1, 3, 7, 16, 33, 67, 100, offset_n }) |n| {
        const m = w * n;
        for (0..8) |ox| for (0..8) |oy| {
            const x = mem + ox;
            const y = mem + w * (offset_n + 8) + oy;
            @memcpy(x[0..m], @as([*]const F, @ptrCast(&vals[0]))[0..m]);
            @memcpy(y[0..m], @as([*]const F, @ptrCast(&vals[1]))[0..m]);
            const ni: BlasLong = @intCast(n);
            _ = switch (op) {
                .scal => @call(.auto, kern, .{ ni, 0, 0 } ++ util.scalarArgs(T, a) ++ .{ x, 1, x, 0, x, 1 }),
                .axpy => @call(.auto, kern, .{ ni, 0, 0 } ++ util.scalarArgs(T, a) ++ .{ x, 1, y, 1, y, 0 }),
                .axpby => @call(.auto, kern, .{ni} ++ util.scalarArgs(T, a) ++ .{ x, 1 } ++ util.scalarArgs(T, b) ++ .{ y, 1 }),
                .rot => kern(ni, x, 1, y, 1, @as(F, 0.6), @as(F, 0.8)),
            };
            if (ox == 0 and oy == 0) {
                @memcpy(first[0][0..m], x[0..m]);
                @memcpy(first[1][0..m], y[0..m]);
            } else if (!std.mem.eql(u8, std.mem.sliceAsBytes(first[0][0..m]), std.mem.sliceAsBytes(x[0..m])) or
                !std.mem.eql(u8, std.mem.sliceAsBytes(first[1][0..m]), std.mem.sliceAsBytes(y[0..m])))
            {
                check.say("{s} {s} n={d} alpha {any}: shifted by {d},{d} gave other bits than unshifted: ", .{ routine, @typeName(T), n, a, ox, oy });
                return error.TestUnexpectedResult;
            }
        };
    }
}

pub fn checkShifted(vt: *const Vtable) !void {
    ed.failed = null;
    inline for (vtable_types) |p| {
        const l = p[0];
        const T = p[1];
        const prec: portable.Prec = comptime .from(l);
        inline for (comptime std.enums.values(Reduce)) |r| check.keep(&ed.failed, shifted(T, portable.slotName(prec, r.op()), false, portable.slotPtr(vt, prec, r.op()).*.?));
        inline for (if (comptime isComplex(T)) .{ portable.Op.dotu_k, portable.Op.dotc_k } else .{portable.Op.dot_k}) |r| check.keep(&ed.failed, shifted(T, portable.slotName(prec, r), true, portable.slotPtr(vt, prec, r).*.?));
        const a: T = if (comptime isComplex(T)) .init(1.5, -0.75) else 1.5;
        // complex kernels branch on a zero part (`scalReal` passes one)
        const scals: []const T = if (comptime isComplex(T)) &.{ a, .init(1.5, 0), .init(0, -0.75) } else &.{a};
        for (scals) |s| check.keep(&ed.failed, shiftedWise(T, l ++ "scal", .scal, portable.slotPtr(vt, .from(l), .scal_k).*.?, s));
        inline for (if (comptime isComplex(T)) .{ portable.Op.axpy_k, portable.Op.axpyc_k } else .{portable.Op.axpy_k}) |r| check.keep(&ed.failed, shiftedWise(T, portable.slotName(prec, r), .axpy, portable.slotPtr(vt, prec, r).*.?, a));
        if (@hasField(Vtable, portable.slotName(.from(l), .axpby_k))) if (portable.slotPtr(vt, .from(l), .axpby_k).*) |k| check.keep(&ed.failed, shiftedWise(T, l ++ "axpby", .axpby, k, a));
        check.keep(&ed.failed, shiftedWise(T, l ++ "rot", .rot, portable.slotPtr(vt, .from(l), .rot_k).*.?, a));
    }
    check.keep(&ed.failed, shifted(f32, "dsdot", true, vt.dsdot_k.?));
    if (ed.failed) |e| return e;
}
