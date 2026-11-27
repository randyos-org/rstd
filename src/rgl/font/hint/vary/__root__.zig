//! Variable TrueType fonts at an instance, in FreeType 2.14.3's fixed
//! point (`src/truetype/ttgxvar.c`, issue #809): `avar` normalization,
//! `gvar` with its inferred deltas, `cvar` and `HVAR`. `MVAR` and `VVAR`
//! are not applied; no test font has them.

pub const Blend = @import("Blend.zig");
pub const Gvar = @import("Gvar.zig");
pub const Hvar = @import("Hvar.zig");
pub const ItemStore = @import("ItemStore.zig");
pub const cvar = @import("cvar.zig");
pub const iup = @import("iup.zig");
pub const packed_data = @import("packed.zig");
pub const tuples = @import("tuples.zig");

test {
    _ = @import("_test/vary.zig");
}
