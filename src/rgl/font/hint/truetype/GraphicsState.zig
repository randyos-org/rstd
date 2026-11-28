//! The graphics state that `fpgm` and `prep` set and every glyph program
//! starts from: FreeType 2.14.3 `ttobjs.h` `TT_GraphicsState` and its
//! defaults, `ttinterp.c` `tt_default_graphics_state`.

const round = @import("round.zig");

/// Reference points.
rp0: u16 = 0,
rp1: u16 = 0,
rp2: u16 = 0,
/// The zone numbers `SZP0` to `SZP2` last set; a run starts with its zone
/// pointers on the glyph zone, these left as they were.
gep0: u16 = 1,
gep1: u16 = 1,
gep2: u16 = 1,
/// 2.14 unit vectors.
dual: [2]i32 = .{ 0x4000, 0 },
proj: [2]i32 = .{ 0x4000, 0 },
free: [2]i32 = .{ 0x4000, 0 },
loop: i32 = 1,
/// The mode last set; a run starts rounding to the grid whatever this says.
round_state: round.Mode = .grid,
/// 26.6 distances.
min_distance: i32 = 64,
control_value_cutin: i32 = 68,
single_width_cutin: i32 = 0,
single_width_value: i32 = 0,
delta_base: u16 = 9,
delta_shift: u16 = 3,
auto_flip: bool = true,
instruct_control: u8 = 0,
scan_control: bool = false,
scan_type: i32 = 0,
