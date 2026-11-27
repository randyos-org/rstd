/* hintrefs <font>: FreeType's outlines of <font>'s test glyphs at the test
 * sizes, unhinted and in each hinting mode, as `src/rgl/font/hint/_test/
 * refs.zig` reads them, on stdout (rstd issue #805). Built and run by
 * tools/nix/hintrefs.sh at a pinned nixpkgs; `zig build hint-refs
 * -Dtest-fonts` regenerates the checked-in files. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_DRIVER_H
#include FT_FONT_FORMATS_H
#include FT_MODULE_H
#include FT_TRUETYPE_TABLES_H

/* ASCII, then Latin-1 (composite glyphs in TrueType fonts), Greek, Cyrillic, the euro sign; Nerd Font icons
 * (no script covers them: the autohinter's CJK fallback style, issue #814); kana, kanji and CJK punctuation;
 * Limbu (the Indic writing system, issue #815) */
static const uint32_t extra[] = {0xC0,   0xC5,   0xC7,   0xC9,   0xDF,   0xE7,   0xE9,   0xF1,   0xF6,   0x3A9,  0x3BB,
                                 0x416,  0x44F,  0x20AC, 0xE0A0, 0xE0B0, 0xF015, 0xF07B, 0xF121, 0x3001, 0x3002, 0x3042,
                                 0x306E, 0x30A2, 0x30AB, 0x4E00, 0x4E2D, 0x56FD, 0x5B57, 0x65E5, 0x6587, 0x672C, 0x6C38,
                                 0x6771, 0x8A9E, 0xFF01, 0x1901, 0x1904, 0x1909, 0x190E, 0x1915, 0x191C, 0x1923, 0x1930,
                                 0x1946, 0x194B};
static const unsigned sizes[] = {9, 12, 16, 24};
#define N_SIZES (sizeof sizes / sizeof sizes[0])

/* refs.zig's `Mode`, by value */
enum { unhinted, native, native_v35, native_light, autohint, autohint_light, n_modes };
static const FT_Int32 mode_flags[n_modes] = {
    FT_LOAD_NO_HINTING,
    FT_LOAD_NO_AUTOHINT,
    FT_LOAD_NO_AUTOHINT,
    FT_LOAD_NO_AUTOHINT | FT_LOAD_TARGET_LIGHT,
    FT_LOAD_FORCE_AUTOHINT,
    FT_LOAD_FORCE_AUTOHINT | FT_LOAD_TARGET_LIGHT,
};
#define no_mode 255

typedef struct {
    FT_Error err;
    FT_Pos advance, lsb_delta, rsb_delta;
    unsigned n_contours, n_points;
    unsigned short *ends;
    unsigned char *tags;
    FT_Vector *points;
} Glyph;

typedef struct {
    unsigned char *p;
    size_t len, cap;
} Buf;

static void die(const char *what, FT_Error err) {
    fprintf(stderr, "hintrefs: %s failed (FreeType error %d)\n", what, err);
    exit(1);
}

static void put(Buf *b, unsigned char c) {
    if (b->len == b->cap) {
        b->cap = b->cap ? 2 * b->cap : 1 << 16;
        b->p = realloc(b->p, b->cap);
        if (!b->p) die("realloc", 0);
    }
    b->p[b->len++] = c;
}

static void varint(Buf *b, uint64_t v) {
    while (v >= 0x80) {
        put(b, (unsigned char)(v | 0x80));
        v >>= 7;
    }
    put(b, (unsigned char)v);
}

static void zigzag(Buf *b, int64_t v) { varint(b, ((uint64_t)v << 1) ^ (uint64_t)(v >> 63)); }

static FT_Face open_face(FT_Library lib, const char *path, unsigned interpreter) {
    FT_Face face;
    FT_Error err = FT_Property_Set(lib, "truetype", "interpreter-version", &interpreter);
    if (err) die("setting the interpreter version", err);
    if ((err = FT_New_Face(lib, path, 0, &face))) die("opening the font", err);
    return face;
}

static Glyph load(FT_Face face, FT_UInt gid, FT_Int32 flags) {
    Glyph g = {0};
    g.err = FT_Load_Glyph(face, gid, flags | FT_LOAD_NO_BITMAP);
    if (g.err) return g;
    FT_GlyphSlot s = face->glyph;
    FT_Outline *o = &s->outline;
    g.advance = s->advance.x;
    g.lsb_delta = s->lsb_delta;
    g.rsb_delta = s->rsb_delta;
    g.n_contours = (unsigned)o->n_contours;
    g.n_points = (unsigned)o->n_points;
    g.ends = malloc(g.n_contours * sizeof *g.ends + 1);
    g.tags = malloc(g.n_points + 1);
    g.points = malloc(g.n_points * sizeof *g.points + 1);
    if (!g.ends || !g.tags || !g.points) die("malloc", 0);
    for (unsigned i = 0; i < g.n_contours; i++) g.ends[i] = (unsigned short)o->contours[i];
    for (unsigned i = 0; i < g.n_points; i++) g.tags[i] = (unsigned char)(FT_CURVE_TAG(o->tags[i]));
    memcpy(g.points, o->points, g.n_points * sizeof *g.points);
    return g;
}

static int same(const Glyph *a, const Glyph *b) {
    if (a->err || b->err) return a->err == b->err;
    return a->advance == b->advance && a->lsb_delta == b->lsb_delta && a->rsb_delta == b->rsb_delta &&
           a->n_contours == b->n_contours && a->n_points == b->n_points &&
           !memcmp(a->ends, b->ends, a->n_contours * sizeof *a->ends) && !memcmp(a->tags, b->tags, a->n_points) &&
           !memcmp(a->points, b->points, a->n_points * sizeof *a->points);
}

static void write_glyph(Buf *b, const Glyph *g) {
    zigzag(b, g->err);
    if (g->err) return;
    zigzag(b, g->advance);
    zigzag(b, g->lsb_delta);
    zigzag(b, g->rsb_delta);
    varint(b, g->n_contours);
    for (unsigned i = 0; i < g->n_contours; i++) varint(b, i ? g->ends[i] - g->ends[i - 1] : g->ends[0]);
    varint(b, g->n_points);
    for (unsigned i = 0; i < g->n_points; i++) put(b, g->tags[i]);
    FT_Pos x = 0, y = 0;
    for (unsigned i = 0; i < g->n_points; i++) {
        zigzag(b, g->points[i].x - x);
        zigzag(b, g->points[i].y - y);
        x = g->points[i].x;
        y = g->points[i].y;
    }
}

static FT_ULong table_len(FT_Face face, const char *tag) {
    FT_ULong len = 0;
    FT_Load_Sfnt_Table(face, FT_MAKE_TAG(tag[0], tag[1], tag[2], tag[3]), 0, NULL, &len);
    return len;
}

/* the mode FreeType's default load (`flags` a target only) matches on every glyph and size, else no_mode */
static unsigned default_mode(FT_Library lib, const char *path, Glyph *ref, unsigned n_glyphs, const FT_UInt *gids,
                             FT_Int32 flags, unsigned a, unsigned b) {
    FT_Face face = open_face(lib, path, 40);
    int match_a = 1, match_b = 1;
    for (unsigned s = 0; s < N_SIZES; s++) {
        FT_Error err = FT_Set_Pixel_Sizes(face, 0, sizes[s]);
        if (err) die("setting the size", err);
        for (unsigned i = 0; i < n_glyphs; i++) {
            if (!gids[i]) continue;
            Glyph g = load(face, gids[i], flags);
            match_a &= same(&g, &ref[(i * N_SIZES + s) * n_modes + a]);
            match_b &= same(&g, &ref[(i * N_SIZES + s) * n_modes + b]);
            free(g.ends), free(g.tags), free(g.points);
        }
    }
    FT_Done_Face(face);
    return match_a ? a : match_b ? b : no_mode;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: hintrefs <font>\n");
        return 2;
    }
    const char *path = argv[1];
    FT_Library lib;
    FT_Error err = FT_Init_FreeType(&lib);
    if (err) die("FT_Init_FreeType", err);

    uint32_t cps[160];
    unsigned n_glyphs = 0;
    for (uint32_t c = 0x20; c < 0x7F; c++) cps[n_glyphs++] = c;
    for (unsigned i = 0; i < sizeof extra / sizeof extra[0]; i++) cps[n_glyphs++] = extra[i];

    FT_Face face = open_face(lib, path, 40);
    FT_UInt gids[160];
    for (unsigned i = 0; i < n_glyphs; i++) gids[i] = FT_Get_Char_Index(face, cps[i]);
    int is_cff = !strcmp(FT_Get_Font_Format(face), "CFF");
    FT_UShort upem = face->units_per_EM;
    FT_ULong fpgm = table_len(face, "fpgm"), prep = table_len(face, "prep"), cvt = table_len(face, "cvt ");
    FT_Done_Face(face);

    Glyph *ref = calloc((size_t)n_glyphs * N_SIZES * n_modes, sizeof *ref);
    if (!ref) die("calloc", 0);
    for (unsigned m = 0; m < n_modes; m++) {
        if (m == native_v35 && is_cff) continue;
        face = open_face(lib, path, m == native_v35 ? 35 : 40);
        for (unsigned s = 0; s < N_SIZES; s++) {
            if ((err = FT_Set_Pixel_Sizes(face, 0, sizes[s]))) die("setting the size", err);
            for (unsigned i = 0; i < n_glyphs; i++)
                if (gids[i]) ref[(i * N_SIZES + s) * n_modes + m] = load(face, gids[i], mode_flags[m]);
        }
        FT_Done_Face(face);
    }
    unsigned def = default_mode(lib, path, ref, n_glyphs, gids, FT_LOAD_DEFAULT, native, autohint);
    unsigned def_light = default_mode(lib, path, ref, n_glyphs, gids, FT_LOAD_TARGET_LIGHT, native_light, autohint_light);

    Buf b = {0};
    const char magic[] = "RSTDHNT1";
    for (unsigned i = 0; i < 8; i++) put(&b, (unsigned char)magic[i]);
    FT_Int major, minor, patch;
    FT_Library_Version(lib, &major, &minor, &patch);
    varint(&b, (uint64_t)major), varint(&b, (uint64_t)minor), varint(&b, (uint64_t)patch);
    varint(&b, upem);
    put(&b, (unsigned char)is_cff);
    varint(&b, fpgm), varint(&b, prep), varint(&b, cvt);
    put(&b, (unsigned char)def), put(&b, (unsigned char)def_light);
    varint(&b, N_SIZES);
    for (unsigned s = 0; s < N_SIZES; s++) varint(&b, sizes[s]);
    varint(&b, n_glyphs);
    for (unsigned i = 0; i < n_glyphs; i++) {
        varint(&b, cps[i]);
        varint(&b, gids[i]);
        if (!gids[i]) continue;
        for (unsigned s = 0; s < N_SIZES; s++)
            for (unsigned m = 0; m < n_modes; m++)
                if (!(m == native_v35 && is_cff)) write_glyph(&b, &ref[(i * N_SIZES + s) * n_modes + m]);
    }
    if (fwrite(b.p, 1, b.len, stdout) != b.len) die("writing", 0);
    FT_Done_FreeType(lib);
    return 0;
}
