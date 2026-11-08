//! Probes `@Matrix` (the randy-matrix fork) through the SPIR-V backend, one fragment entry per use: `zig build-obj <this> -target
//! spirv32-vulkan -mcpu vulkan_v1_2 -fno-llvm -ofmt=spirv -femit-bin=matrix.spv`, then `spirv-val --target-env vulkan1.2 matrix.spv`.
//! `M` is not square and its 3-row columns pad to 16 bytes, so a rows/cols mix-up or a stride bug changes a type.

const std = @import("std");

const M = @Matrix(3, 2, f32);
const Column = @Vector(3, f32);
/// `M`'s array form, `[cols][rows]f32`.
const Columns = [2][3]f32;

comptime {
    const info = @typeInfo(M).matrix;
    std.debug.assert(info.rows == 3 and info.cols == 2 and info.child == f32);
    std.debug.assert(@sizeOf(M) == 2 * @sizeOf(Column));
    std.debug.assert(@alignOf(M) == @alignOf(Column));
}

/// Matrices in a push-constant block, where the column stride is
/// decorated and checked.
const Globals = extern struct {
    xform: M,
    bias: Column,
    /// An array of matrices takes its stride decorations on this member.
    views: [2]@Matrix(4, 4, f32),
};

const out_color = @extern(*addrspace(.output) @Vector(4, f32), .{
    .name = "frag_color",
    .decoration = .{ .location = 0 },
});
const in_a = @extern(*addrspace(.input) Column, .{
    .name = "in_a",
    .decoration = .{ .location = 0 },
});
const in_b = @extern(*addrspace(.input) Column, .{
    .name = "in_b",
    .decoration = .{ .location = 1 },
});
/// Takes locations 2 and 3, one per column.
const in_m = @extern(*addrspace(.input) M, .{
    .name = "in_m",
    .decoration = .{ .location = 2 },
});
const globals = @extern(*addrspace(.push_constant) const Globals, .{
    .name = "globals",
});

/// Writes the first and last row of both columns, so either axis being
/// swapped or truncated shows up in the disassembly.
fn store(cols: Columns) void {
    out_color.* = .{ cols[0][0], cols[0][2], cols[1][0], cols[1][2] };
}

/// A matrix as a function parameter.
fn storeMatrix(m: M) void {
    store(m);
}

/// A matrix as a return value.
fn swapColumns(m: M) M {
    const cols: Columns = m;
    return @as(Columns, .{ cols[1], cols[0] });
}

/// Runtime [cols][rows]T to matrix and back.
fn fromArray() callconv(.{ .spirv_fragment = .{} }) void {
    const a: [3]f32 = in_a.*;
    const b: [3]f32 = in_b.*;
    const arr: Columns = .{ a, b };
    const m: M = arr;
    store(m);
}

/// Runtime tuple of columns to matrix.
fn fromTuple() callconv(.{ .spirv_fragment = .{} }) void {
    const columns = .{ @as([3]f32, in_a.*), @as([3]f32, in_b.*) };
    const m: M = columns;
    store(m);
}

/// A matrix as a shader input (extern matrices are SPIR-V only), loaded
/// into a function-local variable.
fn externInput() callconv(.{ .spirv_fragment = .{} }) void {
    var m: M = undefined;
    m = in_m.*;
    store(m);
}

/// A comptime-known matrix reaching codegen as a constant argument.
fn constant() callconv(.{ .spirv_fragment = .{} }) void {
    const identity: M = @as(Columns, .{ .{ 1, 0, 0 }, .{ 0, 1, 0 } });
    storeMatrix(identity);
}

/// A matrix passed through a function and returned.
fn roundTrip() callconv(.{ .spirv_fragment = .{} }) void {
    storeMatrix(swapColumns(in_m.*));
}

/// A matrix read from a push-constant block.
fn pushConstant() callconv(.{ .spirv_fragment = .{} }) void {
    const cols: Columns = globals.xform;
    const bias: [3]f32 = globals.bias;
    const view = globals.views[1];
    out_color.* = .{ cols[0][0] + bias[0], cols[0][2], cols[1][0], cols[1][2] + bias[2] + view[3][3] };
}

/// Init syntax from runtime columns, then a runtime-indexed column read.
fn initAndIndex() callconv(.{ .spirv_fragment = .{} }) void {
    const m: M = .{ in_a.*, in_b.* };
    var c: usize = @intFromFloat(in_a.*[0]);
    c &= 1;
    const col = m[c];
    out_color.* = .{ col[0], col[2], m[1][0], @floatFromInt(m.len) };
}

/// Element and column stores through a matrix pointer, then a `for` over
/// the columns and a destructure.
fn storeAndLoop() callconv(.{ .spirv_fragment = .{} }) void {
    var m: M = in_m.*;
    m[0][1] = in_a.*[0];
    m[1] = in_b.*;
    var sum: Column = @splat(0);
    for (m) |col| sum += col;
    const c0, _ = m;
    out_color.* = .{ sum[0], sum[1], sum[2], c0[1] };
}

/// Every product form lowers to SPIR-V's own instruction: OpMatrixTimesVector,
/// OpVectorTimesMatrix, OpMatrixTimesMatrix and OpMatrixTimesScalar.
fn products() callconv(.{ .spirv_fragment = .{} }) void {
    const view = globals.views[0];
    const m: M = in_m.*;
    const p: @Vector(4, f32) = .{ in_a.*[0], in_a.*[1], in_a.*[2], 1 };
    const mv = view * p;
    const vm = in_b.* * m;
    const mm = view * globals.views[1];
    const ms = m * in_a.*[0];
    out_color.* = mv + mm[3] + @Vector(4, f32){ vm[0], vm[1], ms[1][2], 0 };
}

/// OpTranspose, and `+`, `-` and negation column by column.
fn transposeAndSums() callconv(.{ .spirv_fragment = .{} }) void {
    const m: M = in_m.*;
    const t: @Matrix(2, 3, f32) = @transpose(m);
    const s = m + m - -m;
    out_color.* = .{ t[2][1], s[1][2], 0, 1 };
}

/// A product outside SPIR-V's matrix sizes falls back to column arithmetic.
fn wideProduct() callconv(.{ .spirv_fragment = .{} }) void {
    const w: @Matrix(4, 5, f32) = .{ globals.views[0][0], globals.views[0][1], globals.views[0][2], globals.views[0][3], globals.views[1][0] };
    const v: @Vector(5, f32) = @splat(in_a.*[0]);
    out_color.* = w * v;
}

comptime {
    @export(&products, .{ .name = "matrix_products" });
    @export(&transposeAndSums, .{ .name = "matrix_transpose_and_sums" });
    @export(&wideProduct, .{ .name = "matrix_wide_product" });
    @export(&initAndIndex, .{ .name = "matrix_init_and_index" });
    @export(&storeAndLoop, .{ .name = "matrix_store_and_loop" });
    @export(&fromArray, .{ .name = "matrix_from_array" });
    @export(&fromTuple, .{ .name = "matrix_from_tuple" });
    @export(&externInput, .{ .name = "matrix_extern_input" });
    @export(&constant, .{ .name = "matrix_constant" });
    @export(&roundTrip, .{ .name = "matrix_round_trip" });
    @export(&pushConstant, .{ .name = "matrix_push_constant" });
}
