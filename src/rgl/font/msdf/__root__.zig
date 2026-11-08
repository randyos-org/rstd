//! Multichannel signed distance fields from glyph outlines, ported from mist
//! (and msdfgen through it; see ../sfnt/readme.md).

pub const Generator = @import("Generator.zig");
pub const error_correction = @import("error_correction.zig");
const zigstd = @import("../../../zigstd.zig");

pub const FontMetrics = Generator.FontMetrics;
pub const GlyphMetricsEm = Generator.GlyphMetricsEm;
pub const KerningPair = Generator.KerningPair;
pub const GeneratedGlyph = Generator.GeneratedGlyph;
pub const GeneratedAtlasGlyph = Generator.GeneratedAtlasGlyph;
pub const GeneratedAtlas = Generator.GeneratedAtlas;
pub const SdfType = Generator.SdfType;
pub const ColoringMethod = Generator.ColoringMethod;
pub const Winding = Generator.Winding;
pub const VarFontArgument = Generator.VarFontArgument;
pub const Options = Generator.Options;
pub const Msdf10Pixel = Generator.Msdf10Pixel;

test {
    zigstd.testRefAllDecls(@This());
}
