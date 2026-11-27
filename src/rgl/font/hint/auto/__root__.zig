//! The autohinter: FreeType 2.14.3's `src/autofit/`, after fontations
//! `skrifa/src/outline/autohint/` (issue #373).

pub const GlyphHints = @import("GlyphHints.zig");
pub const segments = @import("segments.zig");
pub const edges = @import("edges.zig");

test {
    _ = @import("_test/topo.zig");
}
