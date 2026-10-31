//! Probes `@Matrix` (zig-build's zig-src randy-matrix, phase 1) through
//! the SPIR-V backend: one fragment entry point per supported use.
//!
//!   zig build-obj tools/spirv_matrix_probe.zig -target spirv32-vulkan \
//!       -mcpu vulkan_v1_2 -fno-llvm -ofmt=spirv -femit-bin=matrix.spv
//!   spirv-val --target-env vulkan1.2 matrix.spv
//!
//! Phase 1 lowers a matrix as an array of column vectors, not
//! OpTypeMatrix, and has no `.{ col, col }` init syntax, indexing or
//! operators; those are later phases, so this probe should grow with them.
//!
//! `M` is deliberately not square, so a rows/cols mix-up anywhere changes
//! a type, and its 3-row columns are padded to 16 bytes, so the column
//! stride differs from the column's element count.

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
    out_color.* = .{ cols[0][0] + bias[0], cols[0][2], cols[1][0], cols[1][2] + bias[2] };
}

comptime {
    @export(&fromArray, .{ .name = "matrix_from_array" });
    @export(&fromTuple, .{ .name = "matrix_from_tuple" });
    @export(&externInput, .{ .name = "matrix_extern_input" });
    @export(&constant, .{ .name = "matrix_constant" });
    @export(&roundTrip, .{ .name = "matrix_round_trip" });
    @export(&pushConstant, .{ .name = "matrix_push_constant" });
}
