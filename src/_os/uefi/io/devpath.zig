//! Multi-instance device path lists, as the `ConOut` and `ErrOut` variables
//! hold them (UEFI 2.10 sections 3.3 and 10.3.1): instances of nodes, each
//! ended by an end-instance (0x7F, 0x01) or end-entire (0x7F, 0xFF) node.

const std = @import("std");

const end_type = 0x7F;
const end_instance = 0x01;
const end_entire = 0xFF;
const header_len = 4;

const Instances = struct {
    list: []const u8,
    pos: usize = 0,
    done: bool = false,

    /// The next instance, without its end node; null after the end-entire
    /// node, or at a malformed node (`done` stays false then).
    fn next(it: *Instances) ?[]const u8 {
        if (it.done) return null;
        const start = it.pos;
        var at = start;
        while (at + header_len <= it.list.len) {
            const len = std.mem.readInt(u16, it.list[at + 2 ..][0..2], .little);
            if (len < header_len or at + len > it.list.len) return null;
            if (it.list[at] == end_type and it.list[at + 1] == end_instance) {
                it.pos = at + len;
                return it.list[start..at];
            }
            if (it.list[at] == end_type and it.list[at + 1] == end_entire) {
                it.done = true;
                return if (at == start) null else it.list[start..at];
            }
            at += len;
        }
        return null;
    }
};

/// Every node is well formed up to an end-entire node.
fn wellFormed(list: []const u8) bool {
    var it: Instances = .{ .list = list };
    while (it.next()) |_| {}
    return it.done;
}

/// Each instance of `a` is one of `b`'s, and `a` has one.
fn contains(b: []const u8, a: []const u8) bool {
    var ia: Instances = .{ .list = a };
    var any = false;
    while (ia.next()) |x| {
        any = true;
        var ib: Instances = .{ .list = b };
        while (ib.next()) |y| {
            if (std.mem.eql(u8, x, y)) break;
        } else return false;
    }
    return any;
}

/// Each device `b` names is one `a` names, in any order; false if either
/// is empty or malformed.
pub fn covers(a: []const u8, b: []const u8) bool {
    return wellFormed(a) and wellFormed(b) and contains(a, b);
}
