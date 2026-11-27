//! The CVT or storage area a glyph program sees: the size's own until it
//! writes, then a copy, so no glyph's writes reach the next (FreeType 2.14.3
//! `ttinterp.c` `Modify_CVT_Check` and `Ins_WS`'s `glyfStorage`; fontations
//! `skrifa/src/outline/glyf/hint/cow_slice.rs`, commit 4f115482, MIT).

const CowSlice = @This();

/// The size's values, which `fpgm` and `prep` write directly.
base: []i32,
/// Room for a glyph's copy, as long as `base`.
copy: []i32,
/// Whether the glyph program has written, so `copy` is the one in use.
copied: bool = false,

pub fn len(s: CowSlice) usize {
    return s.base.len;
}

pub fn get(s: CowSlice, i: usize) i32 {
    return if (s.copied) s.copy[i] else s.base[i];
}

/// Sets `i`; in a glyph program the first write copies `base` first.
pub fn set(s: *CowSlice, glyph: bool, i: usize, v: i32) void {
    if (glyph and !s.copied) {
        @memcpy(s.copy[0..s.base.len], s.base);
        s.copied = true;
    }
    if (s.copied) s.copy[i] = v else s.base[i] = v;
}

/// Sets `i` in whichever is in use, never copying: FreeType's `Ins_WCVTF`,
/// which skips `Modify_CVT_Check`, so a glyph's WCVTF before any other
/// write reaches the size's CVT.
pub fn setInPlace(s: *CowSlice, i: usize, v: i32) void {
    if (s.copied) s.copy[i] = v else s.base[i] = v;
}
