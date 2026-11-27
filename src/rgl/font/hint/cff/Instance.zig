//! A CFF font at one size: FreeType 2.14.3's scale for it
//! (`FT_Request_Metrics`, `metrics.x_scale` and `y_scale`).

const fixed = @import("../fixed.zig");
const Sfnt = @import("../../sfnt/Sfnt.zig");
const CffTable = @import("../../sfnt/tables/cff/CffTable.zig");

const Instance = @This();

pub const Error = error{ NotCff, InvalidPpem };

sfnt: *const Sfnt,
cff: CffTable,
/// Font units to 26.6 pixels, 16.16.
scale: i32,

pub fn init(sfnt: *const Sfnt, ppem: u16) Error!Instance {
    const cff = sfnt.cff orelse return error.NotCff;
    const head = sfnt.head orelse return error.NotCff;
    if (ppem == 0) return error.InvalidPpem;
    return .{ .sfnt = sfnt, .cff = cff, .scale = fixed.div(@as(i32, ppem) << 6, head.unitsPerEm()) };
}
