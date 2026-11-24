//! The caller's quality-of-service class on macOS, which a worker raises
//! itself to and keeps (issue #691): macOS places a thread by its class, and a
//! default-class worker of a user-interactive caller lands on an Efficiency
//! core. Never above the caller's class; nothing elsewhere.

const sysinfo = @import("builtin");

const Qos = @This();

const on = sysinfo.os.tag == .macos;

extern "c" fn qos_class_self() c_uint;
extern "c" fn pthread_set_qos_class_self_np(qos_class: c_uint, relative_priority: c_int) c_int;

/// `QOS_CLASS_UNSPECIFIED`: nothing to raise to.
const unspecified: c_uint = 0;

class: c_uint = unspecified,

/// The calling thread's class.
pub fn capture() Qos {
    if (comptime !on) return .{};
    return .{ .class = qos_class_self() };
}

/// Raises the calling thread to `q`'s class when that is higher, and keeps it:
/// one task is too short for macOS to move a thread, so a class given back
/// after each does nothing (issue #691). A thread whose class cannot change (a
/// Dispatch worker) keeps its own.
pub fn raise(q: Qos) void {
    if (comptime !on) return;
    if (q.class == unspecified or qos_class_self() >= q.class) return;
    _ = pthread_set_qos_class_self_np(q.class, 0);
}
