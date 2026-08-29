const std = @import("std");
const renderdoc = @import("renderdoc");
const dll = @import("../dll.zig");

pub const RenderDocApi = renderdoc.RENDERDOC_API_1_7_0;
const pRENDERDOC_GetAPI = renderdoc.pRENDERDOC_GetAPI;

pub fn loadRenderDocAPI() ?*RenderDocApi {
    const lib_name = if (@import("builtin").os.tag == .windows) "renderdoc.dll" else "librenderdoc.so";
    // Check if RenderDoc is already injected into this process
    var rdc_lib = dll.findLoaded(lib_name) orelse return null;
    const get_api_fn = rdc_lib.lookup(pRENDERDOC_GetAPI, "RENDERDOC_GetAPI") orelse return null;
    var api_pointers: ?*anyopaque = null;
    if (get_api_fn.?(renderdoc.eRENDERDOC_API_Version_1_7_0, &api_pointers) == 1) {
        return @ptrCast(@alignCast(api_pointers));
    }
    return null;
}

// pub fn main() !void {
//     if (loadRenderDocAPI()) |rdc| {
//         std.debug.print("RenderDoc detected! Triggering a programmatic capture.\n", .{});

//         rdc.StartFrameCapture(null, null);

//         // --- YOUR ZIG RENDER CODE HERE ---
//         // e.g., drawTriangle();

//         _ = rdc.EndFrameCapture(null, null);
//     } else {
//         std.debug.print("Running normally without RenderDoc.\n", .{});
//     }
// }
