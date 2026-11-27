//! A TrueType font at one size, ready to hint its glyphs: FreeType 2.14.3
//! `src/truetype/ttobjs.c` (`tt_size_reset`, `tt_size_init_bytecode`,
//! `tt_size_run_fpgm`, `tt_size_run_prep`, `TT_Save_Context`) and
//! `ttgload.c` `tt_loader_init`: `fpgm` once, `prep` on the scaled CVT,
//! then what every glyph program starts from.

const std = @import("std");
const fixed = @import("../fixed.zig");
const Sfnt = @import("../../sfnt/Sfnt.zig");
const Engine = @import("Engine.zig");
const GraphicsState = @import("GraphicsState.zig");
const definitions = @import("definitions.zig");
const zone = @import("zone.zig");

const Allocator = std.mem.Allocator;
const Instance = @This();

pub const Options = struct {
    ppem: u16,
    version: Engine.Version = .v40,
    target: Engine.Target = .normal,
    pedantic: bool = false,
};

pub const InitError = error{ NotTrueType, InvalidPpem } || Allocator.Error;

allocator: Allocator,
sfnt: *const Sfnt,
options: Options,
/// Font units to 26.6 pixels, 16.16 (`metrics.x_scale` and `y_scale`).
scale: i32,
engine: Engine,
/// What each glyph program starts from (`size->GS`): the defaults and
/// the settings `prep` left.
gs: GraphicsState = .{},
/// `fpgm`'s or `prep`'s error, which fails every hinted load.
failed: ?Engine.Error = null,
/// `prep` turned hinting off (INSTCTRL bit 1).
hinting_off: bool = false,
/// v40's backward compatibility for each load (`tt_loader_init`).
backward_compatibility: u8 = 0,
/// `hdmx`'s widths at this size, in whole pixels, when FreeType uses them.
widths: ?[]const u8 = null,
buffers: []i32,
defs: []definitions.Definition,
twilight: []zone.Point,
twilight_tags: []u8,

/// Sets up `sfnt` at `options.ppem` and runs `fpgm` and `prep`; their
/// errors are kept in `failed`, not returned, as FreeType keeps them.
pub fn init(allocator: Allocator, sfnt: *const Sfnt, options: Options) InitError!Instance {
    const maxp = sfnt.maxp orelse return error.NotTrueType;
    const head = sfnt.head orelse return error.NotTrueType;
    if (sfnt.glyf == null or sfnt.loca == null) return error.NotTrueType;
    if (options.ppem == 0) return error.InvalidPpem;
    const upem = head.unitsPerEm();
    const scale = fixed.div(@as(i32, options.ppem) << 6, upem);
    const max_stack: usize = maxp.maxStackElements();
    const stack_len = max_stack + @max(max_stack / 2, 128);
    const store_len: usize = maxp.maxStorage();
    const cvt_raw = sfnt.cvt orelse &.{};
    const cvt_len = cvt_raw.len / 2;
    const buffers = try allocator.alloc(i32, stack_len + 2 * store_len + 2 * cvt_len);
    errdefer allocator.free(buffers);
    const nf: usize = maxp.maxFunctionDefs();
    const ni: usize = maxp.maxInstructionDefs();
    const defs = try allocator.alloc(definitions.Definition, nf + ni);
    errdefer allocator.free(defs);
    const n_twilight = @as(usize, maxp.maxTwilightPoints()) + 4;
    const twilight = try allocator.alloc(zone.Point, 2 * n_twilight);
    errdefer allocator.free(twilight);
    const twilight_tags = try allocator.alloc(u8, n_twilight);
    @memset(twilight, .{ 0, 0 });
    @memset(twilight_tags, 0);
    var at: usize = 0;
    const take = struct {
        fn f(b: []i32, pos: *usize, n: usize) []i32 {
            defer pos.* += n;
            return b[pos.*..][0..n];
        }
    }.f;
    const stack = take(buffers, &at, stack_len);
    const store = take(buffers, &at, store_len);
    const store_copy = take(buffers, &at, store_len);
    const cvt = take(buffers, &at, cvt_len);
    const cvt_copy = take(buffers, &at, cvt_len);
    var inst: Instance = .{
        .allocator = allocator,
        .sfnt = sfnt,
        .options = options,
        .scale = scale,
        .engine = .{
            .version = options.version,
            .target = options.target,
            .pedantic = options.pedantic,
            // at 72 dpi the point size is the ppem (`tt_size_request`)
            .metrics = .{ .ppem = options.ppem, .scale = scale, .point_size = @as(i32, options.ppem) * 64 },
            .glyph_count = maxp.numGlyphs(),
            .stack = stack,
            .functions = .{ .defs = defs[0..nf] },
            .instructions = .{ .defs = defs[nf..] },
            .cvt = .{ .base = cvt, .copy = cvt_copy },
            .storage = .{ .base = store, .copy = store_copy },
            .zones = .{ .{ .org = twilight[0..n_twilight], .cur = twilight[n_twilight..], .tags = twilight_tags }, .{} },
            .code = .{ sfnt.fpgm orelse &.{}, sfnt.prep orelse &.{}, &.{} },
        },
        .buffers = buffers,
        .defs = defs,
        .twilight = twilight,
        .twilight_tags = twilight_tags,
    };
    inst.runFontProgram();
    if (inst.failed == null) inst.runControlValueProgram(cvt_raw);
    inst.finish();
    return inst;
}

/// An instance whose unhinted loads are in font units (`FT_LOAD_NO_SCALE`);
/// its hinted loads mean nothing.
pub fn initUnscaled(allocator: Allocator, sfnt: *const Sfnt) InitError!Instance {
    var inst = try init(allocator, sfnt, .{ .ppem = 1 });
    inst.scale = 1 << 16;
    inst.hinting_off = true;
    return inst;
}

pub fn deinit(inst: *Instance) void {
    inst.allocator.free(inst.twilight_tags);
    inst.allocator.free(inst.twilight);
    inst.allocator.free(inst.defs);
    inst.allocator.free(inst.buffers);
}

/// `tt_size_run_fpgm`: definitions only; what it leaves of the graphics
/// state `prep` resets.
fn runFontProgram(inst: *Instance) void {
    const e = &inst.engine;
    e.gs = .{};
    e.zones[1] = .{};
    e.run(.font) catch |err| {
        inst.failed = err;
    };
}

/// `tt_size_run_prep`: defaults, zeroed twilight points and storage, the
/// CVT scaled from font units (`FT_MulFix`), then `prep`; its graphics
/// state saved as `TT_Save_Context` saves it.
fn runControlValueProgram(inst: *Instance, cvt_raw: []const u8) void {
    const e = &inst.engine;
    @memset(inst.twilight, .{ 0, 0 });
    @memset(e.storage.base, 0);
    for (e.cvt.base, 0..) |*v, i| {
        const fword: i16 = @bitCast(std.mem.readInt(u16, cvt_raw[2 * i ..][0..2], .big));
        v.* = fixed.mul(fword, inst.scale);
    }
    e.gs = .{};
    e.zones[1] = .{};
    e.run(.control_value) catch |err| {
        inst.failed = err;
        return;
    };
    const g = e.gs;
    inst.gs = .{
        .min_distance = g.min_distance,
        .control_value_cutin = g.control_value_cutin,
        .single_width_cutin = g.single_width_cutin,
        .single_width_value = g.single_width_value,
        .delta_base = g.delta_base,
        .delta_shift = g.delta_shift,
        .auto_flip = g.auto_flip,
        .instruct_control = g.instruct_control,
        .scan_control = g.scan_control,
        .scan_type = g.scan_type,
    };
}

/// `tt_loader_init` after `prep`: INSTCTRL bit 1 turns hinting off, bit 2
/// restores the default graphics state; v40 sets backward compatibility
/// unless the target is mono (tricky fonts too, which this does not
/// detect); `hdmx` gives the advances without it, outside fixed pitch.
fn finish(inst: *Instance) void {
    inst.hinting_off = inst.gs.instruct_control & 1 != 0;
    if (inst.gs.instruct_control & 2 != 0) inst.gs = .{};
    if (inst.options.version == .v40 and inst.options.target != .mono)
        inst.backward_compatibility = (inst.gs.instruct_control & 4) ^ 4;
    const fixed_pitch = if (inst.sfnt.post) |p| p.isFixedPitch() != 0 else false;
    if (inst.backward_compatibility == 0 and !fixed_pitch) inst.widths = hdmxWidths(inst.sfnt, inst.options.ppem);
}

/// `hdmx`'s record for `ppem` (`tt_face_get_device_metrics`), if any.
fn hdmxWidths(sfnt: *const Sfnt, ppem: u16) ?[]const u8 {
    const t = sfnt.hdmx orelse return null;
    const glyphs: usize = (sfnt.maxp orelse return null).numGlyphs();
    if (t.len < 8 or std.mem.readInt(u16, t[0..2], .big) != 0) return null;
    const count = std.mem.readInt(i16, t[2..4], .big);
    const size = std.mem.readInt(u32, t[4..8], .big);
    if (count < 0 or count > 255 or size < glyphs + 2 or size > 0x40000) return null;
    const n: usize = @intCast(count);
    const sz: usize = size;
    if (t.len < 8 + n * sz) return null;
    for (0..n) |r| {
        const rec = t[8 + r * sz ..][0..sz];
        if (rec[0] == ppem) return rec[2..][0..glyphs];
    }
    return null;
}
