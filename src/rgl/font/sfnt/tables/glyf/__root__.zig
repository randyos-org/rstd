//! The 'glyf' table: TrueType outlines. Not decoded up front (a font has
//! thousands of glyphs, callers need a few per frame): `GlyfTable.getGlyph`
//! decodes one glyph's stream on demand from the byte range `LocaTable`
//! gives, and `shape.glyphShapeTT` walks composites into renderable vertices.

pub const GlyfTable = @import("GlyfTable.zig");
pub const SimpleGlyph = @import("SimpleGlyph.zig");
pub const CompositeGlyph = @import("CompositeGlyph.zig");
pub const shape = @import("shape.zig");
