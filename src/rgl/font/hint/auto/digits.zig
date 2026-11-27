//! Whether a font's digits share one advance: FreeType 2.14.3
//! `src/autofit/aflatin.c` (`af_latin_metrics_check_digits`) and `afcjk.c`
//! (`af_cjk_metrics_check_digits`, the same code).

const std = @import("std");
const Source = @import("Source.zig");
const shaper = @import("shaper.zig");

/// True if every digit the font has has the first one's advance.
pub fn sameWidth(src: *Source, allocator: std.mem.Allocator) bool {
    const digits = "0 1 2 3 4 5 6 7 8 9";
    var at: usize = 0;
    var started = false;
    var old_advance: i32 = 0;
    // a failed `FT_Get_Advance` leaves the last value
    var advance: i32 = 0;
    while (at < digits.len) {
        const c = shaper.cluster(src.sfnt, digits, &at);
        if (c.count > 1 or c.glyph == 0) continue;
        if (src.advance(allocator, c.glyph)) |a| advance = a;
        if (started) {
            if (advance != old_advance) return false;
        } else {
            old_advance = advance;
            started = true;
        }
    }
    return true;
}
