//! A font's hinting at any size: none, its designer hints (TrueType
//! bytecode, CFF stem hints), the autohinter, or FreeType 2.14.3's default
//! choice between them (`src/base/ftobjs.c` `FT_Load_Glyph`): the designer
//! hints, unless the load is light and the face TrueType, or the TrueType
//! face has no `fpgm` and a `prep` of at most 7 bytes; never the autohinter
//! for a tricky face. Each size's TrueType instance (`fpgm` and `prep`
//! run) is made once and kept. A variable TrueType font is hinted at the
//! instance `setDesignCoordinates` sets (FreeType's
//! `FT_Set_Var_Design_Coordinates`), else at its default.

const std = @import("std");
const fixed = @import("fixed.zig");
const Sfnt = @import("../sfnt/Sfnt.zig");
const Glyph = @import("Glyph.zig");
const truetype = @import("truetype/__root__.zig");
const cff = @import("cff/__root__.zig");
const auto = @import("auto/__root__.zig");
const tricky_fonts = @import("tricky.zig");
const vary = @import("vary/__root__.zig");

const Allocator = std.mem.Allocator;
const Font = @This();

pub const Mode = enum {
    /// Outlines only scaled (`FT_LOAD_NO_HINTING`).
    none,
    /// The font's own hints (`FT_LOAD_NO_AUTOHINT`).
    designer,
    /// The autohinter (`FT_LOAD_FORCE_AUTOHINT`).
    auto,
    /// FreeType's choice (`FT_LOAD_DEFAULT`).
    default,
};

/// The rendering target, as `FT_LOAD_TARGET_*`.
pub const Target = enum { normal, light };

pub const Options = struct {
    mode: Mode = .default,
    target: Target = .normal,
};

pub const Error = vary.Blend.Error || truetype.Loader.Error || truetype.Instance.InitError || cff.LoadError || cff.Instance.Error || auto.Hinter.Error || auto.Source.InitError;

allocator: Allocator,
sfnt: *const Sfnt,
tricky: bool,
autohinter: ?auto.Hinter = null,
/// The variable font's instance, null at the default.
blend: ?*vary.Blend = null,
/// TrueType instances by size and target.
instances: std.AutoHashMapUnmanaged(Key, *truetype.Instance) = .empty,

const Key = struct { ppem: u16, target: Target };

pub fn init(allocator: Allocator, sfnt: *const Sfnt) Font {
    return .{ .allocator = allocator, .sfnt = sfnt, .tricky = tricky_fonts.isTricky(sfnt) };
}

pub fn deinit(f: *Font) void {
    f.reset();
    if (f.blend) |b| {
        b.deinit();
        f.allocator.destroy(b);
    }
}

/// Drops what was made for the instance.
fn reset(f: *Font) void {
    var it = f.instances.valueIterator();
    while (it.next()) |inst| {
        inst.*.deinit();
        f.allocator.destroy(inst.*);
    }
    f.instances.deinit(f.allocator);
    f.instances = .empty;
    if (f.autohinter) |*a| a.deinit();
    f.autohinter = null;
}

/// Hints at the instance of `design` (16.16 axis values in `fvar` order,
/// missing ones at their defaults); `error.NotVariable` without `fvar`.
pub fn setDesignCoordinates(f: *Font, design: []const i32) vary.Blend.Error!void {
    const next = try vary.Blend.init(f.allocator, f.sfnt, design);
    f.reset();
    if (f.blend) |b| {
        b.deinit();
        f.allocator.destroy(b);
        f.blend = null;
    }
    var n = next orelse return;
    errdefer n.deinit();
    const b = try f.allocator.create(vary.Blend);
    b.* = n;
    f.blend = b;
}

/// What `FT_Load_Glyph` does with a default load of this face.
pub fn defaultMode(f: *const Font, target: Target) Mode {
    if (f.tricky) return .designer;
    // the CFF engine hints lightly itself; the TrueType one does not
    if (f.sfnt.cff == null) {
        if (target == .light) return .auto;
        const fpgm = if (f.sfnt.fpgm) |t| t.len else 0;
        const prep = if (f.sfnt.prep) |t| t.len else 0;
        if (f.sfnt.glyf != null and fpgm == 0 and prep <= 7) return .auto;
    }
    return .designer;
}

fn instance(f: *Font, ppem: u16, target: Target) Error!*truetype.Instance {
    const key: Key = .{ .ppem = ppem, .target = target };
    if (f.instances.get(key)) |i| return i;
    const inst = try f.allocator.create(truetype.Instance);
    errdefer f.allocator.destroy(inst);
    inst.* = try .init(f.allocator, f.sfnt, .{ .ppem = ppem, .target = if (target == .light) .light else .normal, .blend = f.blend });
    errdefer inst.deinit();
    try f.instances.put(f.allocator, key, inst);
    return inst;
}

/// Glyph `gid` at `ppem` in 26.6 pixels; the caller owns it.
pub fn load(f: *Font, allocator: Allocator, gid: u16, ppem: u16, opts: Options) Error!Glyph {
    var mode = opts.mode;
    // a tricky face is never autohinted
    if (mode == .default or (mode == .auto and f.tricky)) mode = f.defaultMode(opts.target);
    switch (mode) {
        .none, .designer => {
            const hinted = mode == .designer;
            if (f.sfnt.cff != null) {
                const inst: cff.Instance = try .init(f.sfnt, ppem);
                return cff.load(allocator, &inst, @fromBackingInt(gid), hinted);
            }
            return truetype.Loader.load(allocator, try f.instance(ppem, opts.target), gid, hinted);
        },
        .auto => {
            if (f.autohinter == null) f.autohinter = try .init(f.allocator, f.sfnt, f.blend);
            return f.autohinter.?.load(allocator, gid, ppem, if (opts.target == .light) .light else .normal);
        },
        .default => unreachable,
    }
}
