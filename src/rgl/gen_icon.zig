const std = @import("std");

const rstd = @import("rstd");
const RGBA = rstd.rgl.color.RGBA;
const RgbaImage = rstd.rgl.RgbaImage;
const downsample = rstd.rgl.draw.downsample;
const codecs = rstd.rgl.codecs;
const png = codecs.png;
const ico = codecs.ico;

pub fn main(init: std.process.Init) !void {
    const alloc = init.arena.allocator();
    const io = init.io;

    const args = try init.minimal.args.toSlice(alloc);
    if (args.len != 5) {
        std.debug.print("usage: genicon <in.png> <out.png> <out.ico> <out.rgba>\n", .{});
        return error.InvalidArgs;
    }
    const in_png_path = args[1];
    const png_path = args[2];
    const ico_path = args[3];
    const raw_path = args[4];

    const png_data = try std.Io.Dir.cwd().readFileAlloc(io, in_png_path, alloc, .unlimited);
    const m = try png.decode(alloc, png_data, .{});

    // Emit PNG (256) and a multi-size ICO.
    const icon256 = try downsample(alloc, &m, 256);
    const png_bytes = try png.encode(alloc, &icon256);
    try std.Io.Dir.cwd().writeFile(io, .{ .sub_path = png_path, .data = png_bytes });

    const sizes = [_]u32{ 16, 24, 32, 48, 64, 128, 256 };
    var imgs: [sizes.len]RgbaImage = undefined;
    for (sizes, 0..) |s, i| imgs[i] = try downsample(alloc, &m, s);

    try ico.writeICO(alloc, io, ico_path, &imgs);

    const raw32 = for (sizes, imgs) |s, im| {
        if (s == 32) break im;
    } else unreachable;
    try raw32.writeRawFile(io, raw_path);
}
