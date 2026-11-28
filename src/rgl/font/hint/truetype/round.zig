//! The interpreter's rounding modes (FreeType 2.14.3 `ttinterp.c`, the
//! `Round_*` functions, their compensation always 0).

const fixed = @import("../fixed.zig");

pub const Mode = enum { half_grid, grid, double_grid, down_to_grid, up_to_grid, off, super, super45 };

/// `SROUND` and `S45ROUND`'s period, phase and threshold, in 26.6.
pub const Super = struct {
    period: i32 = 64,
    phase: i32 = 0,
    threshold: i32 = 0,
};

/// `d` rounded as `mode` rounds it, keeping its sign.
pub fn round(mode: Mode, s: Super, d: i32) i32 {
    const neg = d < 0;
    // the magnitude, rounded; FreeType's `SUB_LONG( compensation, d )`
    const m = if (neg) 0 -% d else d;
    switch (mode) {
        .off => return d,
        .grid => return clamp(neg, fixed.round(m), 0),
        .half_grid => return clamp(neg, fixed.floor(m) +% 32, 32),
        .double_grid => return clamp(neg, fixed.roundPad(m, 32), 0),
        .down_to_grid => return clamp(neg, fixed.floor(m), 0),
        .up_to_grid => return clamp(neg, fixed.ceil(m), 0),
        .super => {
            const bias = s.threshold -% s.phase;
            return clamp(neg, ((m +% bias) & -%s.period) +% s.phase, s.phase);
        },
        .super45 => {
            const bias = s.threshold -% s.phase;
            const q = if (s.period == -1) 0 -% (m +% bias) else @divTrunc(m +% bias, s.period);
            return clamp(neg, (q *% s.period) +% s.phase, s.phase);
        },
    }
}

/// `v` with `d`'s sign: `-v` for negative `d`, and `floor` (signed) when
/// the result crossed zero, as each `Round_*` falls back.
fn clamp(neg: bool, v: i32, floor: i32) i32 {
    if (!neg) return if (v < 0) floor else v;
    const r = 0 -% v;
    return if (r > 0) -%floor else r;
}
