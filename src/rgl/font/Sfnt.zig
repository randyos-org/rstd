const tables = @import("tables.zig");

const Self = @This();

cmap: tables.CmapTable,
head: tables.HeadTable,
hhea: tables.HheaTable,
hmtx: tables.HmtxTable,

/// only required if not CFF
glyf: ?tables.GlyphTable = null,
/// only required if glyf is present
loca: ?tables.LocaTable = null,
/// only required if not glyf
cff: ?tables.CffTable = null,

/// optional
maxp: ?tables.MaxpTable = null,
