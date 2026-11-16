//! Each shaderlib probe (`src/rgl/shaders/probe.zig`) drawn on the device for
//! a few inputs, its one pixel read back as `RGBA32F` and checked against the
//! same code on the CPU (`probe.compute`, issue #380). The driver sets each
//! instruction's precision and may fuse `a*b + c` (Vulkan's precision rules,
//! <https://registry.khronos.org/vulkan/specs/latest/html/vkspec.html#spirvenv-precision-operation>;
//! design.md, "Bits under SPIR-V"), so each op has an error budget in ulps of
//! its inputs' magnitude rather than a bit promise.

const std = @import("std");
const vk = @import("vulkan");
const rstd = @import("rstd");
const probe = @import("probe");
const Headless = @import("Headless.zig");

const gpu = rstd.rgl.vulkan.gpu;
const Allocator = std.mem.Allocator;
const Op = probe.Op;
const Vec4 = @Vector(4, f32);

/// The vertex attributes the probes' vertex stage reads, locations 0 to 8.
const Instance = extern struct {
    a: [3]f32,
    b: [3]f32,
    c: [3]f32,
    t: f32,
    m: [4][4]f32,
    v: [4]f32,

    fn of(in: probe.Inputs) Instance {
        return .{ .a = in.a, .b = in.b, .c = in.c, .t = in.t, .m = .{ in.m[0], in.m[1], in.m[2], in.m[3] }, .v = in.v };
    }
};

/// A copy of SPIR-V `bytes` on a word boundary, as `vkCreateShaderModule`
/// reads words and `@embedFile` promises bytes.
fn aligned(gpa: Allocator, bytes: []const u8) ![]align(4) u8 {
    const copy = try gpa.alignedAlloc(u8, .@"4", bytes.len);
    @memcpy(copy, bytes);
    return copy;
}

fn attr(location: u32, format: vk.Format, offset: usize) vk.VertexInputAttributeDescription {
    return .{ .binding = 0, .location = location, .format = format, .offset = @intCast(offset) };
}

const attrs = [_]vk.VertexInputAttributeDescription{
    attr(0, .r32g32b32_sfloat, @offsetOf(Instance, "a")),
    attr(1, .r32g32b32_sfloat, @offsetOf(Instance, "b")),
    attr(2, .r32g32b32_sfloat, @offsetOf(Instance, "c")),
    attr(3, .r32_sfloat, @offsetOf(Instance, "t")),
    attr(4, .r32g32b32a32_sfloat, @offsetOf(Instance, "m")),
    attr(5, .r32g32b32a32_sfloat, @offsetOf(Instance, "m") + 16),
    attr(6, .r32g32b32a32_sfloat, @offsetOf(Instance, "m") + 32),
    attr(7, .r32g32b32a32_sfloat, @offsetOf(Instance, "m") + 48),
    attr(8, .r32g32b32a32_sfloat, @offsetOf(Instance, "v")),
};

/// Inputs on which every op is defined: `a` positive (`inversesqrt`),
/// `a[0] < b[0]` (`smoothstep`'s edges), `b <= c` (`clamp`'s bounds), `t`
/// positive, `m` diagonally dominant (`inv`).
const cases = [_]probe.Inputs{
    .{
        .a = .{ 0.5, 1.25, 2 },
        .b = .{ 0.75, 1.5, 2.5 },
        .c = .{ 1.5, 1.75, 3 },
        .t = 0.8,
        .m = .{ .{ 4, 1, 0.5, 0.25 }, .{ 0.5, 3, 0.75, 1 }, .{ 0.25, 0.5, 5, 1.5 }, .{ 1, 0.25, 0.5, 2 } },
        .v = .{ 1, -2, 0.5, 3 },
    },
    .{
        .a = .{ 0.25, 3.5, 1.125 },
        .b = .{ 1.75, 0.5, 0.625 },
        .c = .{ 2, 1, 0.875 },
        .t = 1.3,
        .m = .{ .{ -3, 0.5, 1, 0.25 }, .{ 0.75, 4.5, -0.5, 1 }, .{ 0.5, -1, 3.5, 0.75 }, .{ 0.25, 1, -0.5, -2.5 } },
        .v = .{ -0.5, 1.5, 2, -1 },
    },
};

/// `op`'s budget: `k` ulps of the inputs' largest magnitude to the power
/// `p`, the terms an exact result sums; `p` null measures by the result
/// itself (relative error: `inversesqrt`, `inv`).
const Budget = struct { k: f32, p: ?u32 };

fn budget(op: Op) Budget {
    return switch (op) {
        // correctly rounded, or exact
        .step, .clamp, .clamp_vec, .minmax, .fract, .fract_scalar => .{ .k = 2, .p = 1 },
        .mix, .mix_scalar, .faceforward => .{ .k = 4, .p = 1 },
        .smoothstep, .normalize => .{ .k = 16, .p = 0 },
        // products, maybe fused
        .vec_arith, .fma_muladd, .fma_builtin => .{ .k = 4, .p = 2 },
        .dot, .dot_reduce, .dot_asm, .dot2, .cross, .cross_asm, .outer => .{ .k = 8, .p = 2 },
        .mat_vec, .vec_mat => .{ .k = 16, .p = 2 },
        .mat_mat => .{ .k = 64, .p = 2 },
        // through `inversesqrt` or `sqrt`: a few ulps over the products
        .length, .distance => .{ .k = 16, .p = 1 },
        .reflect => .{ .k = 16, .p = 3 },
        .refract => .{ .k = 32, .p = 3 },
        .det => .{ .k = 64, .p = 4 },
        .inversesqrt, .inversesqrt_vec => .{ .k = 4, .p = null },
        .inv => .{ .k = 512, .p = null },
    };
}

fn largest(in: probe.Inputs) f32 {
    var m: f32 = @abs(in.t);
    for ([_][3]f32{ in.a, in.b, in.c }) |x| for (x) |e| {
        m = @max(m, @abs(e));
    };
    for (in.m ++ [_]Vec4{in.v}) |col| m = @max(m, @reduce(.Max, @abs(col)));
    return @max(m, 1);
}

/// `got`'s error over `op`'s budget, for the inputs `in`: within at 1 or
/// less; infinity for a NaN where `want` has none.
fn spent(op: Op, in: probe.Inputs, want: f32, got: f32) f32 {
    if (std.math.isNan(want) or std.math.isNan(got)) return if (std.math.isNan(want) and std.math.isNan(got)) 0 else std.math.inf(f32);
    const b = budget(op);
    const scale = if (b.p) |p| std.math.pow(f32, largest(in), @floatFromInt(p)) else @max(@abs(want), std.math.floatMin(f32));
    return @abs(got - want) / (b.k * std.math.floatEps(f32) * scale);
}

/// Runs every probe on every case; the count of mismatches.
pub fn run(h: *Headless, gpa: Allocator) !u32 {
    const target = try Headless.Target.init(h, .r32g32b32a32_sfloat, 1, 1);
    defer target.deinit(h);
    const out = try h.buffer(@sizeOf(Vec4), .{ .transfer_dst = true });
    defer out.deinit(h);
    var failures: u32 = 0;
    inline for (comptime std.enums.values(Op)) |op| {
        const name = @tagName(op);
        const vert = try aligned(gpa, @embedFile(name ++ ".vert.spv"));
        defer gpa.free(vert);
        const frag = try aligned(gpa, @embedFile(name ++ ".frag.spv"));
        defer gpa.free(frag);
        const pipeline = try gpu.createPipelineWith(h.dev, target.render_pass, vert, frag, &attrs, @sizeOf(Instance), false, .{ .blend = false });
        defer pipeline.deinit(h.dev);
        var worst: f32 = 0;
        for (cases, 0..) |in, ci| {
            const instance = Instance.of(in);
            const buf = try h.bufferWith(std.mem.asBytes(&instance), .{ .vertex_buffer = true });
            defer buf.deinit(h);
            try h.submit(Draw{ .h = h, .target = target, .pipeline = pipeline.handle, .instance = buf.buf, .out = out }, Draw.record);
            const got = std.mem.bytesToValue([4]f32, (try out.map(h))[0..@sizeOf(Vec4)]);
            out.unmap(h);
            const want: [4]f32 = probe.compute(op, in);
            for (0..4) |i| {
                const s = spent(op, in, want[i], got[i]);
                worst = @max(worst, s);
                if (s <= 1) continue;
                failures += 1;
                std.debug.print("probe {s} case {d} [{d}]: got {e}, want {e}\n", .{ name, ci, i, got[i], want[i] });
            }
        }
        std.debug.print("probe {s}: worst error {d:.3} of its budget\n", .{ name, worst });
    }
    return failures;
}

const Draw = struct {
    h: *const Headless,
    target: Headless.Target,
    pipeline: vk.Pipeline,
    instance: vk.Buffer,
    out: Headless.Buffer,

    fn record(d: Draw, cmd: vk.CommandBuffer) !void {
        d.target.begin(d.h, cmd, .{ 0, 0, 0, 0 });
        d.h.dev.cmdBindPipeline(cmd, .graphics, d.pipeline);
        d.h.dev.cmdBindVertexBuffers(cmd, 0, &.{d.instance}, &.{0});
        d.h.dev.cmdDraw(cmd, 3, 1, 0, 0);
        d.target.end(d.h, cmd, d.out);
    }
};
