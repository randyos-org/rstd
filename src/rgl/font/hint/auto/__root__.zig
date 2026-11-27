//! The autohinter: FreeType 2.14.3's `src/autofit/`, after fontations
//! `skrifa/src/outline/autohint/` (issue #373).

pub const GlyphHints = @import("GlyphHints.zig");
pub const segments = @import("segments.zig");
pub const edges = @import("edges.zig");
pub const widths = @import("widths.zig");
pub const Source = @import("Source.zig");
pub const shaper = @import("shaper.zig");
pub const coverage = @import("coverage.zig");
pub const LatinMetrics = @import("latin/Metrics.zig");
pub const scripts = @import("data/scripts.zig");
pub const styles = @import("data/styles.zig");
pub const blues = @import("data/blues.zig");

test {
    _ = @import("_test/topo.zig");
    _ = @import("_test/metrics.zig");
}
