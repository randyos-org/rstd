// Fetched via the calling shader's own root rather than
// `@import("shader_opts")` directly to avoid circular dependency issues.
const shader_opts = @import("root").shader_opts;

const types = @import("../types.zig");
const ShaderStage = types.ShaderStage;

pub inline fn getCurrentShaderStage() ShaderStage {
    return comptime @field(ShaderStage, @tagName(shader_opts.stage));
}

pub const VertFn = fn () callconv(.spirv_vertex) void;
pub const FragFn = fn () callconv(.{ .spirv_fragment = .{} }) void;

fn StageFn() type {
    return switch (getCurrentShaderStage()) {
        .vert => VertFn,
        .frag => FragFn,
    };
}

/// `const main = exportShaderMain(struct { pub const vert = v; pub const frag = f; });`
/// `comptime { _ = main; }`
/// Struct decls are lazy, so only the current stage's fn gets compiled;
/// the other stage's stores to `param`s would be read-only input writes.
pub inline fn exportShaderMain(comptime S: type) *const StageFn() {
    const main: *const StageFn() = switch (getCurrentShaderStage()) {
        .vert => &S.vert,
        .frag => &S.frag,
    };
    @export(main, .{ .name = "main" });
    return main;
}

fn Input(comptime T: type) type {
    return *addrspace(.input) T;
}

fn Output(comptime T: type) type {
    return *addrspace(.output) T;
}

fn Param(comptime T: type) type {
    return switch (getCurrentShaderStage()) {
        .vert => Output(T),
        .frag => Input(T),
    };
}

pub const frag_color_out = @extern(*addrspace(.output) types.Rgba, .{
    .name = "frag_color",
    .decoration = .{ .location = 0 },
});

/// For vertex it is output, for fragment it is input.
pub inline fn param(comptime name: []const u8, comptime T: type, location: u32) Param(T) {
    comptime return @extern(Param(T), .{ .name = name, .decoration = .{ .location = location } });
}

pub inline fn paramFlat(comptime name: []const u8, comptime T: type, location: u32) Param(T) {
    comptime return @extern(Param(T), .{ .name = name, .decoration = .{ .flat = location } });
}

pub inline fn input(comptime name: []const u8, comptime T: type, location: u32) Input(T) {
    comptime return @extern(Input(T), .{ .name = name, .decoration = .{ .location = location } });
}

const ConstSampledImage = *addrspace(.constant) const types.SampledImage;

pub inline fn constImage(comptime name: []const u8, set: u32, binding: u32) ConstSampledImage {
    comptime return @extern(ConstSampledImage, .{
        .name = name,
        .decoration = .{ .descriptor = .{ .set = set, .binding = binding } },
    });
}
