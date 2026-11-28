//! A font's three programs (fontations `read-fonts/src/tables/glyf/bytecode/program.rs`, commit 4f115482, MIT).

/// Which program code comes from: only `font` may define functions and
/// instructions, and only `glyph` has a glyph to hint.
pub const Program = enum(u2) {
    /// `fpgm`
    font,
    /// `prep`
    control_value,
    /// a glyph's own
    glyph,
};
