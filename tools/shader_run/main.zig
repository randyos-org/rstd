//! `zig build shader-run -Dvulkan`: shaderlib's probes and rgl's pipelines
//! run on a headless Vulkan device, the first with a graphics queue (Mesa's
//! lavapipe in CI), and checked (issue #380). Exits 1 on any mismatch.

const std = @import("std");
const Headless = @import("Headless.zig");
const probes = @import("probes.zig");
const draws = @import("draws.zig");

pub fn main(init: std.process.Init) !void {
    const gpa = init.gpa;
    var h: Headless = undefined;
    try h.init(gpa);
    defer h.deinit();
    std.debug.print("shader-run on {s}\n", .{h.deviceName()});
    const probe_failures = try probes.run(&h, gpa);
    std.debug.print("probes: {d} mismatches\n", .{probe_failures});
    const draw_failures = try draws.run(&h, gpa);
    std.debug.print("draws: {d} mismatches\n", .{draw_failures});
    if (probe_failures + draw_failures != 0) std.process.exit(1);
}
