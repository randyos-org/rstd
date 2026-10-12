//! Microsoft Office file formats.

/// [MS-CFB] compound file (OLE structured storage) writer.
pub const cfb = @import("cfb.zig");
/// [MS-OVBA] VBA project (`vbaProject.bin`) writer.
pub const ovba = @import("ovba.zig");
/// SpreadsheetML (`.xlsx`/`.xlsm` part) helpers.
pub const xlsx = @import("xlsx.zig");

test {
    const zigstd = @import("../zigstd.zig");
    zigstd.testRefAllDecls(@This());
}
