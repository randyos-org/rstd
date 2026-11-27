//! Glyph hinting that emulates FreeType 2.14.3 (issue #373): designer hints
//! where a font has them, the autohinter where it has none, in FreeType's
//! fixed-point arithmetic so every target gets the same outlines.

const zigstd = @import("../../../zigstd.zig");

pub const fixed = @import("fixed.zig");
pub const truetype = @import("truetype/__root__.zig");
pub const cff = @import("cff/__root__.zig");
pub const auto = @import("auto/__root__.zig");
pub const trig = @import("trig.zig");
pub const Glyph = @import("Glyph.zig");
pub const Font = @import("Font.zig");
pub const tricky = @import("tricky.zig");

test {
    zigstd.testRefAllDecls(@This());
    _ = @import("_test/refs_check.zig");
    _ = @import("_test/truetype_refs.zig");
    _ = @import("_test/cff_refs.zig");
    _ = @import("_test/auto_refs.zig");
    _ = @import("_test/font.zig");
    _ = @import("_test/tricky.zig");
}
