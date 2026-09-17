pub const general = @import("general.zig");
pub const HeadTable = @import("head.zig").HeadTable;
pub const MaxpTable = @import("maxp.zig").MaxpTable;
pub const HheaTable = @import("hhea.zig").HheaTable;
pub const HmtxTable = @import("hmtx.zig").HmtxTable;
pub const LocaTable = @import("loca.zig").LocaTable;
pub const NameTable = @import("name.zig").NameTable;
pub const PostTable = @import("post.zig").PostTable;
pub const KernTable = @import("kern.zig").KernTable;
pub const KernSubtable = @import("kern.zig").KernSubtable;
pub const Os2Table = @import("os2.zig").Os2Table;
pub const FvarTable = @import("fvar.zig").FvarTable;
pub const AvarTable = @import("avar.zig").AvarTable;
pub const GvarTable = @import("gvar.zig").GvarTable;
pub const CvarTable = @import("cvar.zig").CvarTable;
pub const HvarTable = @import("hvar.zig").HvarTable;
pub const VvarTable = @import("vvar.zig").VvarTable;
pub const MvarTable = @import("mvar.zig").MvarTable;
pub const StatTable = @import("stat.zig").StatTable;

const cff = @import("cff/__root__.zig");
pub const CffTable = cff.CffTable;
pub const CffError = cff.CffError;
pub const CharstringCtx = cff.CharstringCtx;
pub const CharstringError = cff.CharstringError;

const gpos = @import("gpos/__root__.zig");
pub const GposTable = gpos.GposTable;

const glyf = @import("glyf/__root__.zig");
pub const GlyfTable = glyf.GlyfTable;

const cmap = @import("cmap/__root__.zig");
pub const CmapTable = cmap.CmapTable;
pub const CmapSubtable = cmap.CmapSubtable;
pub const CmapEncRecord = cmap.EncRecord;

test {
    const std = @import("std");
    std.testing.refAllDecls(@This());
    std.testing.refAllDecls(cff);
    std.testing.refAllDecls(gpos);
    std.testing.refAllDecls(glyf);
    std.testing.refAllDecls(cmap);
}
