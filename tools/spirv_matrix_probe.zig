//! Probes `@Matrix` (zig-build's zig-src e7950e99f, phase 1) through the
//! SPIR-V backend: one fragment entry point per supported use.
//!
//!   zig build-obj tools/spirv_matrix_probe.zig -target spirv32-vulkan \
//!       -mcpu vulkan_v1_2 -fno-llvm -ofmt=spirv -femit-bin=matrix.spv
//!   spirv-val --target-env vulkan1.2 matrix.spv
//!
//! Phase 1 lowers a matrix as an array of column vectors, not
//! OpTypeMatrix, and has no `.{ col, col }` init syntax; both are later
//! phases, so this probe should grow with them.

const M = @Matrix(2, 2, f32);

const out_color = @extern(*addrspace(.output) @Vector(4, f32), .{
    .name = "frag_color",
    .decoration = .{ .location = 0 },
});
const in_a = @extern(*addrspace(.input) @Vector(2, f32), .{
    .name = "in_a",
    .decoration = .{ .location = 0 },
});
const in_b = @extern(*addrspace(.input) @Vector(2, f32), .{
    .name = "in_b",
    .decoration = .{ .location = 1 },
});
const in_m = @extern(*addrspace(.input) M, .{
    .name = "in_m",
    .decoration = .{ .location = 2 },
});

fn store(cols: [2][2]f32) void {
    out_color.* = .{ cols[0][0], cols[0][1], cols[1][0], cols[1][1] };
}

/// The type alone, read at comptime.
fn layout() callconv(.{ .spirv_fragment = .{} }) void {
    const size: f32 = @floatFromInt(@sizeOf(M));
    out_color.* = .{ size, 0, 0, 1 };
}

/// Runtime [cols][rows]T to matrix and back.
fn fromArray() callconv(.{ .spirv_fragment = .{} }) void {
    const a: [2]f32 = in_a.*;
    const b: [2]f32 = in_b.*;
    const arr: [2][2]f32 = .{ a, b };
    const m: M = arr;
    store(m);
}

/// Runtime tuple of columns to matrix.
fn fromTuple() callconv(.{ .spirv_fragment = .{} }) void {
    const columns = .{ @as([2]f32, in_a.*), @as([2]f32, in_b.*) };
    const m: M = columns;
    store(m);
}

/// A matrix as a shader input (extern matrices are SPIR-V only).
fn externInput() callconv(.{ .spirv_fragment = .{} }) void {
    store(in_m.*);
}

comptime {
    @export(&layout, .{ .name = "matrix_layout" });
    @export(&fromArray, .{ .name = "matrix_from_array" });
    @export(&fromTuple, .{ .name = "matrix_from_tuple" });
    @export(&externInput, .{ .name = "matrix_extern_input" });
}
