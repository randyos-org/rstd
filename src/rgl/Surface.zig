//! Backend-agnostic drawing target. Pick a backend once, at creation time;
//! every render.zig draw call afterward is blind to which one it is.
//!
//! sw draws eagerly, straight into `impl.sw`'s pixel buffer. vulkan instead
//! accumulates instance data in `impl.vulkan`; nothing is submitted to a
//! GPU until a consumer that owns the actual device/swapchain (none exists
//! in rstd yet; see rgl/vulkan/gpu.zig) uploads and records it.

const sw = @import("sw.zig");
const vulkan = @import("vulkan/__root__.zig");
const zigstd = @import("../zigstd.zig");
const Allocator = zigstd.Allocator;

const Self = @This();

pub const Backend = enum { sw, vulkan };

pub const Impl = union(Backend) {
    sw: sw.SwSurface,
    vulkan: vulkan.VkSurface,
};

pub const InitOptions = union(Backend) {
    sw: sw.InitOptions,
    vulkan: vulkan.InitOptions,
};

width: u32,
height: u32,
impl: Impl,

pub fn init(alloc: Allocator, width: u32, height: u32, options: InitOptions) !Self {
    const impl: Impl = switch (options) {
        .sw => |o| .{ .sw = try sw.SwSurface.init(alloc, width, height, o) },
        .vulkan => |o| .{ .vulkan = try vulkan.VkSurface.init(alloc, width, height, o) },
    };
    return .{ .width = width, .height = height, .impl = impl };
}

pub fn deinit(self: *Self, alloc: Allocator) void {
    switch (self.impl) {
        inline else => |*s| s.deinit(alloc),
    }
}

test "init/deinit round-trips for the sw backend" {
    const zigstd_test = @import("../zigstd.zig");
    var surface = try Self.init(zigstd_test.testAllocator, 4, 4, .{ .sw = .{} });
    defer surface.deinit(zigstd_test.testAllocator);
    try zigstd_test.testExpectEqual(@as(u32, 4), surface.width);
}
