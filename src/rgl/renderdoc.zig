//! Programmatic RenderDoc captures. Without rstd's `-Drenderdoc` (or on
//! non-desktop targets) every call is a no-op.

const std = @import("std");
const sysinfo = @import("builtin");
const dll = @import("../dll.zig");

const impl = if (@import("rstd_options").has_renderdoc) Api else NoApi;
pub const initRenderDoc = impl.init;
pub const startCapture = impl.start;
pub const endCapture = impl.end;

const NoApi = struct {
    fn init() void {}
    fn start() !void {}
    fn end() !void {}
};

const Api = struct {
    const renderdoc = @import("renderdoc");
    const RenderDocApi = renderdoc.RENDERDOC_API_1_7_0;
    const pRENDERDOC_GetAPI = renderdoc.pRENDERDOC_GetAPI;
    var rdoc_inst: ?*RenderDocApi = null;
    var rdoc_cap_started: bool = false;
    var rdoc_load_failed: bool = false;

    fn loadRenderDocAPI() ?*RenderDocApi {
        const lib_name = if (sysinfo.target.os.tag == .windows) "renderdoc.dll" else "librenderdoc.so";
        // Check if RenderDoc is already injected into this process
        var rdc_lib = dll.findLoaded(lib_name) orelse return null;
        const get_api_fn = rdc_lib.lookup(pRENDERDOC_GetAPI, "RENDERDOC_GetAPI") orelse return null;
        var api_pointers: ?*anyopaque = null;
        if (get_api_fn.?(renderdoc.eRENDERDOC_API_Version_1_7_0, &api_pointers) == 1) {
            return @ptrCast(@alignCast(api_pointers));
        }
        return null;
    }

    fn init() void {
        if (rdoc_inst == null and !rdoc_load_failed) {
            rdoc_inst = loadRenderDocAPI();
            if (rdoc_inst == null) {
                rdoc_load_failed = true;
                std.debug.print("RenderDoc API not available.\n", .{});
            } else {
                std.debug.print("RenderDoc API initialized!\n", .{});
            }
        }
    }

    fn start() !void {
        init();
        if (rdoc_inst) |inst| {
            if (rdoc_cap_started) {
                std.debug.print("RenderDoc capture already started.\n", .{});
                return error.CaptureAlreadyStarted;
            }
            std.debug.print("RenderDoc detected! Triggering a programmatic capture.\n", .{});
            rdoc_cap_started = true;
            inst.StartFrameCapture.?(null, null);
        }
    }

    fn end() !void {
        if (rdoc_inst) |inst| {
            if (!rdoc_cap_started) {
                std.debug.print("RenderDoc capture not started.\n", .{});
                return error.CaptureNotStarted;
            }
            rdoc_cap_started = false;
            if (inst.EndFrameCapture.?(null, null) == 0) {
                std.debug.print("RenderDoc capture failed.\n", .{});
                return error.CaptureFailedToEnd;
            }
            std.debug.print("RenderDoc capture ended successfully.\n", .{});
        }
    }
};
