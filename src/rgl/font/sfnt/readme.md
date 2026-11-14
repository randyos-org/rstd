# rfont - Randy's Font Library

A font library formed by the merger of [sphtext from sphimp](https://github.com/sphaerophoria/sphimp/blob/master/src/sphtext), [TrueType](https://codeberg.org/andrewrk/TrueType) (font parsing), [Mist](https://github.com/flut2/mist) (MSDF generation), and [Turbopack](https://github.com/flut2/turbopack) (atlas packing).  The goal is to create an MSDF-based font rendering library with no libc dependencies, purely implemented in Zig.  It parses `.ttf`/`.otf` files and produces MSDF glyph atlases suitable for scalable, GPU-rendered text.

*Note: For full disclosure, code from sphtext is unlicensed, but a previous message from sphaerophoria in his Discord server requesting clarification on another repo stated "do what you want with it."  That statement is probably not legally binding here, but I'm going to gamble that all his unlicensed code is public domain.  I am also heavily making modifications to the code, so it may end up being largely a re-implementation from scratch after all should this become a problem later.

## How it works

Turning a font file into a renderable glyph atlas goes through several stages. A handful of terms come straight from the OpenType/TrueType file format specification, or from the distance-field literature.  They're also the actual names used in the code (table names, tags, etc.), so it's worth knowing what they mean rather than avoiding them. Paths here are relative to `src/rgl/font/`.

1. **Load** (`sfnt/Sfnt.zig`): every `.ttf`/`.otf` file is wrapped in the same container format, called *sfnt*: a short header followed by a directory of *tables*.  Each table is identified by a 4-byte *tag* (`cmap`, `glyf`, `hhea`, `hmtx`, `maxp`, `kern`, `GPOS`, and `"CFF "` [sic] are all such tags). `init()` reads that directory and records where each table lives in the file. A collection (`.ttc`/`.otc`) prefixes several such directories with a `ttcf` header; `numberOfFonts()` counts them and `init()`'s `font_index` option picks one.

   It also checks which of two mutually exclusive ways the font stores its glyph outlines:

   - `glyf`: TrueType curves; `glyf` is the literal 4-byte tag
   - `"CFF "`: Compact Font Format, a PostScript-derived alternative

   That choice decides which of two independent outline-reading backends every later glyph lookup goes through. CFF fonts additionally get their CFF-specific header parsed (`sfnt/tables/cff/CffTable.zig`) to locate the *charstrings* index (see step 3) and subroutine indices.

2. **Character to glyph mapping** (`sfnt/tables/cmap/`): a font doesn't draw "the letter A" directly; it draws *glyphs*, identified by a small integer called a *glyph index* (several different characters can even share one glyph, or one character can be drawn by different glyphs depending on context). The `cmap` table is what maps Unicode codepoints (`'A'`, `'€'`, ...) to glyph indices.  (The `cmap` table actually supports several internal sub-formats for how that mapping is encoded; formats 0, 2, 4, 6, 8, 10, 12, 13 and 14 are read, one file each.) `Sfnt.getGlyphIndex` searches it and returns a `GlyphIndex`.

3. **Outline extraction** (`Sfnt.glyphShape`, dispatching to `sfnt/tables/glyf/` or `sfnt/tables/cff/charstring.zig`): turns a `GlyphIndex` into a flat list of drawing commands (move-to, line-to, quadratic curve-to, cubic curve-to) describing the glyph's outline, in the font's own internal coordinate units (`[]Vertex`).  Everything above this point is shared by both of the rendering paths--they only differ in what they do with the resulting `[]Vertex`. The two outline formats are read in completely different ways:

   - `glyf`: walks the `glyf` table's raw point/flag/coordinate byte arrays directly, handling both simple contours and *compound glyphs* - glyphs built by combining other glyphs, e.g. an accented letter made from a base letter plus a separately-defined accent mark.
   - `"cff "`: runs a bytecode interpreter. CFF fonts don't store outlines as raw point lists; they store *Type2 charstrings*.  This is a compact, stack-based bytecode format whose instructions mean things like "move by this much," "draw a curve through these points," or "call this reusable subroutine" (subroutines let a font share common shapes, like the curve of a serif, across many glyphs instead of repeating them). `runCharstring` is the small virtual machine that executes that bytecode. It runs twice per glyph: once just to count how many vertices the shape needs, then again to actually fill them in, as an allocation-avoidance strategy instead of using a growable list.

4. **Fast path: direct coverage-bitmap rendering** (`render/rasterize.zig`, `render/tessellate.zig`, `render/bounds.zig`): for cases that just need a glyph rendered once at a fixed pixel size (no scaling later), the vertex outline is rasterized straight to an 8-bit coverage bitmap: each pixel's value says how much of that pixel the glyph's shape covers. `flattenCurves` recursively subdivides every curve into short straight line segments until the approximation is imperceptibly close, producing a plain polygon per contour. `rasterize` then fills those polygons with antialiasing: it computes the exact fractional area of each pixel the shape covers, sweeping scanline by scanline using signed trapezoid/triangle area math. That gives smooth edges without the cost of supersampling (rendering at a higher resolution and downscaling). `rasterize.glyphBitmap` ties this together: get the vertices, compute the pixel-space bounding box (`glyphBitmapBox`), allocate a buffer sized to it, rasterize into it, and return a `GlyphBitmap`: width, height, and an offset describing where that bitmap sits relative to the glyph's logical origin point.

5. **Scalable path: MSDF generation** (`msdf/`, `msdf.Generator`): for glyphs that need to render crisply at any size on the GPU (the primary purpose of this library), the vertex outline instead becomes an MSDF (multichannel signed distance field): a small texture where each texel stores, per color channel, the signed distance from that point to the nearest edge of the glyph outline (negative inside the shape, positive outside). A GPU shader can then reconstruct a crisp edge at *any* zoom level with a single texture sample plus a threshold test, because the distance field itself scales smoothly; unlike a plain coverage bitmap, which just gets blurry or blocky when scaled. Multichannel (as opposed to a plain single-channel signed distance field, or SDF) means the distance is encoded across three channels using different subsets of nearby edges, which is what lets sharp corners get reconstructed sharply instead of rounding off - the classic weakness of single-channel SDFs. Producing one requires a few steps of its own:
   - The vertex outline is converted into contours made of *edge segments* (line, quadratic, or cubic segments) - the same shape, just represented the way the distance-field math operates on it directly, rather than flattened into a polygon like the fast path above.
   - Each edge is assigned a color - some combination of red/green/blue - so that later, each output channel only "sees" a subset of the outline's edges. This is what preserves sharp corners: wherever two differently-colored edges meet, the reconstruction step can tell there's a corner there instead of smoothing it away.
   - For every texel in the output, the true signed distance to the nearest edge of each color is computed analytically (not by rasterizing first) and written to that texel's corresponding channel.
   - A correction pass inspects the generated texels for cases where the per-channel reconstruction would actually produce a visible artifact (an incorrect edge or a strange notch) and fixes those texels directly, since the math above is exact but not automatically artifact-free for the full range of ways outlines can be shaped.

6. **Atlas packing** (`render/atlas.zig`, with `src/rgl/pack.zig`'s rectangle packer): generating one MSDF bitmap per glyph on demand would mean binding a different texture for every glyph drawn, which is far too slow for real text rendering. Instead, the MSDF bitmaps for every glyph a renderer will need get packed together into one shared texture, called an atlas: a rectangle-packing algorithm decides where each glyph's small bitmap goes within the larger shared image, and that placement (plus the glyph's size and positioning metrics) is recorded so that a renderer can look up the right region of the atlas for any glyph with a single texture bind.

Two side paths skip outline extraction entirely, since they're just further table lookups:

- `Sfnt.verticalMetrics` and `Sfnt.getGlyphMetrics` read the `hhea` (horizontal header), `hmtx` (horizontal metrics) and `OS/2` tables for each glyph's *advance width* (how far to move the pen after drawing it) and the font's overall ascent/descent/line-gap.
- `Sfnt.glyphKernAdvance` reads `GPOS` (the modern OpenType Glyph Positioning table, `sfnt/tables/gpos/`) or, if that's absent, the older `kern` table (`sfnt/tables/kern.zig`), to get an extra positioning adjustment for specific pairs of glyphs (e.g. tucking a lowercase "o" slightly under a capital "T" so the spacing looks even).

## Synopsis

A single glyph rendered straight to a coverage bitmap (the fast path):

```zig
const font = rstd.rgl.font;
const Sfnt = font.sfnt.Sfnt;
const sfnt = try Sfnt.init(alloc, @embedFile("GoNotoCurrent-Regular.ttf"), .{});

const glyph = sfnt.getGlyphIndex('A') orelse return error.MissingGlyph;
const scale = sfnt.scaleForPixelHeight(20.0);

var pixels: std.ArrayList(u8) = .empty;
defer pixels.deinit(alloc);
const bmp = try font.render.rasterize.glyphBitmap(sfnt, alloc, &pixels, glyph, scale, scale);
// pixels.items[0 .. bmp.width * bmp.height] is the 8-bit coverage bitmap;
// bmp.off_x/off_y is where it sits relative to the glyph origin.
```

An MSDF atlas for a line of text (the scalable path), from the same `sfnt`:

```zig
var atlas = try font.render.atlas.buildTextAtlas(sfnt, alloc, "Hello, world!", 32, 0, 0, .{ .type = .msdf });
defer atlas.deinit(alloc);
// atlas.pixels is one shared texture (kerning-aware, left to right);
// atlas.glyphs[i].{pos,size,atlas_pos} place each glyph on screen and in it.
```

## Features and Limitations

- Codepoint to glyph lookup, kerning, and both direct coverage-bitmap rendering and MSDF/SDF atlas generation.
- Both TrueType (`glyf`) and CFF/Type2 outline formats are supported.
- Font collections (`.ttc`/`.otc`) load one font at a time by index.
- Font shaping and ligatures are not supported: a run of text is drawn glyph-by-glyph via `cmap` lookups and kerning pairs, not through the more elaborate substitution/positioning rules (as used by, say, Arabic script shaping or ligature substitution like "fi" merging into a single glyph) that a dedicated shaping engine like HarfBuzz would apply.
- Untrusted font files are not supported - the parsers here trust the input to be well-formed.

## Upstream baselines

The upstream commit each source was ported from, and the newest upstream commit already reviewed for changes worth pulling in. When syncing, review from "Checked through" forward, then update that column.

| Source | Ported from | Checked through |
| --- | --- | --- |
| [TrueType](https://codeberg.org/andrewrk/TrueType) | `8e5b3216906ce32463a0db75e3889b2b4df00bf1` (2026-08-20) | `4fea1e7c0f9d4d79348d0302c7628f485301b6ac` (2026-10-29) |
| [mist](https://github.com/flut2/mist) | `eb7b593d685befadc666e3d74cbe8ded6894a3cd` (2026-07-28) | `d44f70b0ee744dde72b726f05b93e96a13a795f9` (2026-10-10) |
| [msdfgen](https://github.com/Chlumsky/msdfgen) | via mist | `1c106ed8117893bf943e577f62eb0665fb271e46` (2026-08-29) |
| [Turbopack](https://github.com/flut2/turbopack) | `90cc3da6bd6f4d5da142651b861a1c6834f3af9c` (2026-05-24) | `90cc3da6bd6f4d5da142651b861a1c6834f3af9c` (2026-05-24) |
| [sphimp](https://github.com/sphaerophoria/sphimp) | `b0c1843652638bb77c354c257c1a730b3e23ba60` (2025-04-02) | `b0c1843652638bb77c354c257c1a730b3e23ba60` (2025-04-02) |

TrueType's baseline is exact (rgui's old submodule pin). The others weren't recorded at import time; each is the upstream head when the code landed here (mist 2026-09-15, the rest by 2026-09-14), and mist's matches our `msdf/error_correction.zig` line for line.
