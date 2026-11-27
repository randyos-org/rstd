//! Whether `system_table.std_err` is a console of its own (#740): stderr
//! goes to it and to `con_out` then, so no device loses the log; else to
//! `con_out` alone, written once.

const std = @import("std");
const uefi = std.os.uefi;

const devpath = @import("devpath.zig");

const SimpleTextOutput = uefi.protocol.SimpleTextOutput;

/// Room for each of the `ConOut` and `ErrOut` variables; a longer one
/// counts as different devices.
const max_variable_len = 1024;

/// `std_err` when it reaches a device `con_out` does not, else null.
pub fn find(con_out: ?*SimpleTextOutput) ?*SimpleTextOutput {
    const st = uefi.system_table;
    const err = st.std_err orelse return null;
    const out = con_out orelse return err;
    if (err == out) return null;
    if (st.standard_error_handle != null and st.standard_error_handle == st.console_out_handle) return null;
    return if (outCoversErr(st.runtime_services)) null else err;
}

/// Each device in the `ErrOut` variable is one in `ConOut`, so `con_out`
/// alone reaches them all.
fn outCoversErr(rs: *uefi.tables.RuntimeServices) bool {
    var out_buf: [max_variable_len]u8 = undefined;
    var err_buf: [max_variable_len]u8 = undefined;
    const out = variable(rs, std.unicode.utf8ToUtf16LeStringLiteral("ConOut"), &out_buf) orelse return false;
    const err = variable(rs, std.unicode.utf8ToUtf16LeStringLiteral("ErrOut"), &err_buf) orelse return false;
    return devpath.covers(out, err);
}

fn variable(rs: *uefi.tables.RuntimeServices, name: [*:0]const u16, buf: []u8) ?[]const u8 {
    const got = rs.getVariable(name, &uefi.tables.global_variable, buf) catch return null;
    return (got orelse return null)[0];
}
