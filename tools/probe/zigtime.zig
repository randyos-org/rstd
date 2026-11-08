//! `zig run tools/probe/zigtime.zig -- [--top N] <zig arguments>`: one compile's `--time-report` (phases, and the declarations with the
//! most `sema_ns`), over the compiler's `--listen=-` protocol. `zig` must be on `PATH`.

const std = @import("std");
const Io = std.Io;
const Server = std.zig.Server;

var A: std.mem.Allocator = undefined;

fn die(comptime f: []const u8, args: anytype) noreturn {
    std.debug.print(f ++ "\n", args);
    std.process.exit(1);
}

/// Reads exactly `buf.len` bytes; false at the end of the stream.
fn readExact(file: Io.File, io: Io, buf: []u8) !bool {
    var got: usize = 0;
    while (got < buf.len) {
        const n = try file.readStreaming(io, &.{buf[got..]});
        if (n == 0) return false;
        got += n;
    }
    return true;
}

const Decl = struct { name: []const u8, file: []const u8, count: u32, sema_ns: u64, codegen_ns: u64, link_ns: u64 };

fn ms(ns: u64) f64 {
    return @as(f64, @floatFromInt(ns)) / 1e6;
}

fn bySema(_: void, a: Decl, b: Decl) bool {
    return a.sema_ns > b.sema_ns;
}

fn report(body: []const u8, top: usize) void {
    const H = Server.Message.TimeReport;
    var h: H = undefined;
    @memcpy(std.mem.asBytes(&h), body[0..@sizeOf(H)]);
    const s = h.stats;
    std.debug.print(
        \\
        \\compiler time (cpu ms): parse {d:.0}, astgen {d:.0}, sema {d:.0}, codegen {d:.0}, link {d:.0}
        \\real ms: files {d:.0}, decls {d:.0}, llvm emit {d:.0}, link flush {d:.0}
        \\reachable files {d}, imported {d}, generic instances {d}, inline calls {d}
        \\
    , .{
        ms(s.cpu_ns_parse),  ms(s.cpu_ns_astgen),   ms(s.cpu_ns_sema),       ms(s.cpu_ns_codegen),     ms(s.cpu_ns_link),
        ms(s.real_ns_files), ms(s.real_ns_decls),   ms(s.real_ns_llvm_emit), ms(s.real_ns_link_flush), s.n_reachable_files,
        s.n_imported_files,  s.n_generic_instances, s.n_inline_calls,
    });
    var at: usize = @sizeOf(H) + h.llvm_pass_timings_len;
    var files: std.ArrayList([]const u8) = .empty;
    for (0..h.files_len) |_| {
        const end = std.mem.findScalarPos(u8, body, at, 0).?;
        files.append(A, body[at..end]) catch @panic("OOM");
        at = end + 1;
    }
    var decls: std.ArrayList(Decl) = .empty;
    for (0..h.decls_len) |_| {
        const end = std.mem.findScalarPos(u8, body, at, 0).?;
        const name = body[at..end];
        at = end + 1;
        const file = std.mem.readInt(u32, body[at..][0..4], .little);
        const count = std.mem.readInt(u32, body[at + 4 ..][0..4], .little);
        const sema = std.mem.readInt(u64, body[at + 8 ..][0..8], .little);
        const cg = std.mem.readInt(u64, body[at + 16 ..][0..8], .little);
        const ln = std.mem.readInt(u64, body[at + 24 ..][0..8], .little);
        at += 32;
        decls.append(A, .{ .name = name, .file = files.items[file], .count = count, .sema_ns = sema, .codegen_ns = cg, .link_ns = ln }) catch @panic("OOM");
    }
    std.mem.sort(Decl, decls.items, {}, bySema);
    std.debug.print("{d} declarations; the {d} with the most semantic analysis:\n", .{ decls.items.len, @min(top, decls.items.len) });
    std.debug.print("{s:>10} {s:>7} {s:>10} {s:>10}  {s}\n", .{ "sema ms", "times", "codegen ms", "link ms", "declaration (file)" });
    for (decls.items[0..@min(top, decls.items.len)]) |d| {
        std.debug.print("{d:>10.1} {d:>7} {d:>10.1} {d:>10.1}  {s} ({s})\n", .{ ms(d.sema_ns), d.count, ms(d.codegen_ns), ms(d.link_ns), d.name, d.file });
    }
}

pub fn main(init: std.process.Init) !void {
    A = init.arena.allocator();
    const io = init.io;
    const args = try init.minimal.args.toSlice(A);
    var top: usize = 25;
    var i: usize = 1;
    if (i + 1 < args.len and std.mem.eql(u8, args[i], "--top")) {
        top = try std.fmt.parseInt(usize, args[i + 1], 10);
        i += 2;
    }
    if (i >= args.len) die("usage: zigtime [--top N] <zig arguments>", .{});
    var argv: std.ArrayList([]const u8) = .empty;
    try argv.append(A, "zig");
    try argv.appendSlice(A, args[i..]);
    try argv.appendSlice(A, &.{ "--time-report", "--listen=-" });

    var child = try std.process.spawn(io, .{
        .argv = argv.items,
        .stdin = .pipe,
        .stdout = .pipe,
        .stderr = .inherit,
    });
    // `update`, no body
    try child.stdin.?.writeStreamingAll(io, &[_]u8{ 1, 0, 0, 0, 0, 0, 0, 0 });

    var reported = false;
    while (true) {
        var hdr: [8]u8 = undefined;
        if (!try readExact(child.stdout.?, io, &hdr)) break;
        const tag = std.mem.readInt(u32, hdr[0..4], .little);
        const len = std.mem.readInt(u32, hdr[4..8], .little);
        const body = try A.alloc(u8, len);
        if (!try readExact(child.stdout.?, io, body)) break;
        const t: Server.Message.Tag = @fromBackingInt(@intCast(tag));
        switch (t) {
            .time_report => {
                report(body, top);
                reported = true;
            },
            // an update ends with the emitted file's digest, or with errors
            .emit_digest => break,
            .error_bundle => {
                std.debug.print("the compilation has errors (see the compiler's own output)\n", .{});
                break;
            },
            else => {},
        }
    }
    // `exit`, no body
    try child.stdin.?.writeStreamingAll(io, &[_]u8{ 0, 0, 0, 0, 0, 0, 0, 0 });
    child.stdin.?.close(io);
    child.stdin = null;
    _ = try child.wait(io);
    if (!reported) std.debug.print("no time report received\n", .{});
}
