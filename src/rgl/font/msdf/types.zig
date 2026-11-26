//! The generator's options and what it returns.

const std = @import("std");
const zigstd = @import("../../../zigstd.zig");
const pack = @import("../../pack.zig");
const error_correction = @import("error_correction.zig");
const Scanline = @import("Scanline.zig");
const coords = @import("../sfnt/coords.zig");
const Allocator = zigstd.Allocator;

pub const FontMetrics = struct {
    line_height: f64,
    ascender: f64,
    descender: f64,
    underline_y: f64,
    underline_thickness: f64,
};

/// Em-scaled
pub const GlyphMetricsEm = struct {
    advance: f64,
    bearing_x: f64,
    bearing_y: f64,
    width: u16,
    height: u16,
};

pub const KerningPair = struct {
    codepoint_1: u21,
    codepoint_2: u21,
    x: f64,
    y: f64,
};

pub const GeneratedGlyph = struct {
    metrics: GlyphMetricsEm,
    pixels: []const u8,

    pub fn deinit(self: GeneratedGlyph, allocator: Allocator) void {
        allocator.free(self.pixels);
    }
};

pub const GeneratedAtlasGlyph = struct {
    glyph_data: GlyphMetricsEm,
    codepoint: u21,
    /// This is the glyph's unscaled, unnormalized
    /// bounding box on the atlas, with padding included.
    tex_bounds: pack.Rect,
};

pub const GeneratedAtlas = struct {
    glyphs: []const GeneratedAtlasGlyph,
    kernings: []const KerningPair,
    /// In case of `msdf10` you should reinterpret this as a `Msdf10Pixel` slice,
    /// like with 10-bit ABGR GPU formats (assuming a little endian host).
    pixels: []const u8,

    pub fn deinit(self: GeneratedAtlas, allocator: Allocator) void {
        allocator.free(self.glyphs);
        allocator.free(self.kernings);
        allocator.free(self.pixels);
    }
};

pub const SdfType = enum {
    sdf,
    psdf,
    msdf,
    mtsdf,
    /// Experimental: A packed BGR MSDF where each channel is 10-bit,
    /// with a 2-bit alpha channel that is ignored (set to u2 max).
    ///
    /// Can prove useful in place of MSDFs as native `R8G8B8_X` (and equivalent)
    /// format support is scarce. Additionally, 3-channel images are often padded
    /// to have an alignment of 4 bytes per pixel on a lot of hardware,
    /// which results in the final byte getting wasted on such formats.
    msdf10,

    pub fn numChannels(self: SdfType) u8 {
        return switch (self) {
            .sdf, .psdf => 1,
            .msdf, .msdf10 => 3,
            .mtsdf => 4,
        };
    }

    /// Bytes per output texel: msdf10 packs its three channels into four.
    pub fn bytesPerPixel(self: SdfType) u8 {
        return switch (self) {
            .msdf10 => 4,
            else => self.numChannels(),
        };
    }

    pub fn requiresColoring(self: SdfType) bool {
        return switch (self) {
            .msdf, .msdf10, .mtsdf => true,
            else => false,
        };
    }
};

pub const ColoringMethod = enum {
    simple,
    /// Only for use with ink trap fonts, as the coloring remains correct
    /// after removing the edges required for trapping ink.
    ink_trap,
    /// Performs the coloring based on edge distances.
    /// Somewhat slower than other methods, but it produces a better result most of the time.
    distance,
};

pub const Winding = enum {
    /// Attempts to figure out winding on its own, by checking
    /// the polarity of an OOB point's distance.
    guess,
    positive,
    negative,
};

/// An axis tag ("wght") and its user-scale value.
pub const VarFontArgument = coords.AxisValue;

pub const Options = struct {
    sdf_type: SdfType,
    px_size: u16,
    px_range: u16,
    /// Has no effect if `sdf_type.requiresColoring()` is false.
    coloring_rng_seed: u64 = 0,
    /// The method with which to perform the MSDF 3-coloring.
    /// While the implementations are based on msdfgen, they're (intentionally)
    /// not equivalent, but should resolve corners similarly well.
    ///
    /// Has no effect if `sdf_type.requiresColoring()` is false.
    coloring_method: ColoringMethod = .distance,
    /// The angle which is considered to be a corner, in radians.
    corner_angle_threshold: f64 = 3.0,
    winding: Winding = .guess,
    /// Validates that the given (or generated) shapes' contours form a
    /// closed loop, with each edge connecting to each other properly.
    validate_shape: bool = false,
    normalize_shape: bool = false,
    orient_contours: bool = false,
    /// Requires `orient_contours` to be disabled.
    scanline_fill_rule: ?Scanline.FillRule = null,
    /// Only MSDFs (both their normal and their 10-bit versions) and MTSDFs can be error corrected.
    error_correction_opts: ?error_correction.Options = .{},
    /// Axis values for a variable font; axes left out take their defaults.
    var_font_args: []const VarFontArgument = &.{},
    /// Whether to use async tasks over concurrent ones during atlas generation.
    /// Currently has no effect outside of atlas generation.
    disable_concurrency: bool = false,
};

pub const Msdf10Pixel = packed struct(u32) {
    r: u10 = 0,
    g: u10 = 0,
    b: u10 = 0,
    a: u2 = std.math.maxInt(u2),
};

/// A distance texel: one value per channel.
pub fn DistPixel(comptime ty: SdfType) type {
    return [ty.numChannels()]f64;
}

/// An output texel, as stored in GeneratedGlyph/GeneratedAtlas pixels.
pub fn OutPixel(comptime ty: SdfType) type {
    return [ty.bytesPerPixel()]u8;
}
