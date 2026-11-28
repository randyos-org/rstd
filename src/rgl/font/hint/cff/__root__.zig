//! CFF hinting: FreeType 2.14.3's Adobe engine (`src/psaux/`), after
//! fontations `skrifa/src/outline/cff/` (issue #373).

pub const Instance = @import("Instance.zig");
pub const load = @import("load.zig").load;
pub const LoadError = @import("load.zig").Error;
pub const charstring = @import("charstring.zig");
pub const private = @import("private.zig");
pub const hintmap = @import("hintmap.zig");
pub const Blues = @import("blues.zig").Blues;
pub const Builder = @import("Builder.zig");
pub const GlyphPath = @import("GlyphPath.zig");

test {
    _ = @import("_test/path.zig");
    _ = @import("_test/seac.zig");
}
