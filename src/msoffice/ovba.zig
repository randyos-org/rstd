//! [MS-OVBA] VBA project writer: builds a `vbaProject.bin` from module
//! source text, with no Office installation needed.
//!
//! Only source code is stored. `_VBA_PROJECT` carries version 0xFFFF and no
//! performance cache, which per [MS-OVBA] 2.3.4.1 tells the host there is no
//! valid compiled p-code, so Office compiles the modules from source on load.
//! Module names and source must be ASCII (a subset of every MBCS code page).
const std = @import("std");
const zigstd = @import("../zigstd.zig");
const bytes = @import("../bytes.zig");
const cfb = @import("cfb.zig");

const Allocator = zigstd.Allocator;
const ByteArrayList = bytes.ByteArrayList;

pub const ModuleKind = enum {
    /// A standard module (`Module=` in PROJECT).
    procedural,
    /// The workbook's document module (ThisWorkbook).
    workbook,
    /// A worksheet's document module.
    worksheet,
};

pub const Module = struct {
    name: []const u8,
    kind: ModuleKind,
    /// Module body: no `Attribute` lines (they are generated), LF or CRLF.
    code: []const u8 = "",
};

pub const Reference = struct {
    name: []const u8,
    /// [MS-OVBA] 2.1.1.8 LibidReference. Hosts resolve the type library by
    /// GUID and version through the registry; the path is only a hint.
    libid: []const u8,
};

/// The references Excel gives every new VBA project. Office is needed for
/// `IRibbonControl` in ribbon callbacks.
pub const default_references = [_]Reference{
    .{
        .name = "stdole",
        .libid = "*\\G{00020430-0000-0000-C000-000000000046}#2.0#0#C:\\Windows\\System32\\stdole2.tlb#OLE Automation",
    },
    .{
        .name = "Office",
        .libid = "*\\G{2DF8D04C-5BFA-101B-BDE5-00AA0044DE52}#2.0#0#C:\\Program Files\\Common Files\\Microsoft Shared\\OFFICE16\\MSO.DLL#Microsoft Office 16.0 Object Library",
    },
};

pub const Project = struct {
    name: []const u8 = "VBAProject",
    /// `{GUID}` identifying the project; also seeds the protection-state
    /// encryption key.
    id: []const u8,
    modules: []const Module,
    references: []const Reference = &default_references,
    code_page: u16 = 1252,
    lcid: u32 = 0x0409,
};

const workbook_base = "0{00020819-0000-0000-C000-000000000046}";
const worksheet_base = "0{00020820-0000-0000-C000-000000000046}";

/// Builds the complete `vbaProject.bin` compound file. Caller owns the result.
pub fn build(gpa: Allocator, project: Project) ![]u8 {
    var arena_state: std.heap.ArenaAllocator = .init(gpa);
    defer arena_state.deinit();
    const arena = arena_state.allocator();

    try checkAscii(project.name);
    try checkAscii(project.id);

    const vba_nodes = try arena.alloc(cfb.Node, project.modules.len + 2);
    vba_nodes[0] = .{ .name = "_VBA_PROJECT", .kind = .{ .stream = &vba_project_stream } };
    vba_nodes[1] = .{ .name = "dir", .kind = .{ .stream = try compress(arena, try dirStream(arena, project)) } };
    for (project.modules, vba_nodes[2..]) |m, *node| {
        try checkAscii(m.name);
        const source = try moduleSource(arena, m);
        node.* = .{ .name = m.name, .kind = .{ .stream = try compress(arena, source) } };
    }

    const root = [_]cfb.Node{
        .{ .name = "PROJECT", .kind = .{ .stream = try projectStream(arena, project) } },
        .{ .name = "PROJECTwm", .kind = .{ .stream = try projectWmStream(arena, project) } },
        .{ .name = "VBA", .kind = .{ .storage = vba_nodes } },
    };
    return cfb.write(gpa, &root);
}

/// [MS-OVBA] 2.3.4.1: Reserved1, Version (0xFFFF on write), Reserved2,
/// Reserved3, and an empty performance cache.
const vba_project_stream = [_]u8{ 0xCC, 0x61, 0xFF, 0xFF, 0x00, 0x03, 0x00 };

fn checkAscii(s: []const u8) !void {
    for (s) |c| if (c >= 0x80) return error.VbaNonAscii;
}

/// Strips a VBE export header (`VERSION`, `BEGIN`..`END`, `Attribute` lines)
/// from the start of exported module text, leaving only the code.
pub fn stripExportHeader(text: []const u8) []const u8 {
    var rest = text;
    var in_begin = false;
    while (rest.len > 0) {
        const eol = std.mem.findScalar(u8, rest, '\n') orelse rest.len;
        const line = std.mem.trimEnd(u8, rest[0..eol], "\r");
        const header = if (in_begin)
            true
        else
            std.mem.startsWith(u8, line, "VERSION ") or
                std.mem.startsWith(u8, line, "Attribute ") or
                zigstd.equal(u8, line, "BEGIN");
        if (!header) break;
        if (zigstd.equal(u8, line, "BEGIN")) in_begin = true;
        if (in_begin and zigstd.equal(u8, line, "END")) in_begin = false;
        rest = if (eol < rest.len) rest[eol + 1 ..] else rest[rest.len..];
    }
    return rest;
}

/// Full module text as stored: generated attributes, then the code, with
/// CRLF line endings throughout.
fn moduleSource(arena: Allocator, m: Module) ![]u8 {
    var out: ByteArrayList = .empty;
    try appendLine(arena, &out, try arena.print("Attribute VB_Name = \"{s}\"", .{m.name}));
    switch (m.kind) {
        .procedural => {},
        .workbook, .worksheet => {
            try appendLine(arena, &out, if (m.kind == .workbook)
                "Attribute VB_Base = \"" ++ workbook_base ++ "\""
            else
                "Attribute VB_Base = \"" ++ worksheet_base ++ "\"");
            for ([_][]const u8{
                "Attribute VB_GlobalNameSpace = False",
                "Attribute VB_Creatable = False",
                "Attribute VB_PredeclaredId = True",
                "Attribute VB_Exposed = True",
                "Attribute VB_TemplateDerived = False",
                "Attribute VB_Customizable = True",
            }) |line| try appendLine(arena, &out, line);
        },
    }
    try checkAscii(m.code);
    var it = std.mem.splitScalar(u8, m.code, '\n');
    while (it.next()) |raw| {
        const line = std.mem.trimEnd(u8, raw, "\r");
        if (it.peek() == null and line.len == 0) break; // trailing newline
        try appendLine(arena, &out, line);
    }
    return out.items;
}

fn appendLine(gpa: Allocator, out: *ByteArrayList, line: []const u8) !void {
    try out.appendSlice(gpa, line);
    try out.appendSlice(gpa, "\r\n");
}

// ---------------------------------------------------------------------------
// PROJECT and PROJECTwm streams ([MS-OVBA] 2.3.1, 2.3.3)
// ---------------------------------------------------------------------------

fn projectStream(arena: Allocator, p: Project) ![]u8 {
    var out: ByteArrayList = .empty;
    try appendLine(arena, &out, try arena.print("ID=\"{s}\"", .{p.id}));
    for (p.modules) |m| {
        const line = switch (m.kind) {
            .procedural => try arena.print("Module={s}", .{m.name}),
            .workbook, .worksheet => try arena.print("Document={s}/&H00000000", .{m.name}),
        };
        try appendLine(arena, &out, line);
    }
    try appendLine(arena, &out, try arena.print("Name=\"{s}\"", .{p.name}));
    try appendLine(arena, &out, "HelpContextID=\"0\"");
    try appendLine(arena, &out, "VersionCompatible32=\"393222000\"");

    // Protection state: not locked, no password, visible ([MS-OVBA] 2.3.1.15-17).
    const key = projectKey(p.id);
    try appendEncryptedLine(arena, &out, "CMG", 0x36, key, &.{ 0, 0, 0, 0 });
    try appendEncryptedLine(arena, &out, "DPB", 0x6C, key, &.{0});
    try appendEncryptedLine(arena, &out, "GC", 0xA2, key, &.{0xFF});

    try appendLine(arena, &out, "");
    try appendLine(arena, &out, "[Host Extender Info]");
    try appendLine(arena, &out, "&H00000001={3832D640-CF90-11CF-8E43-00A0C911005A};VBE;&H00000000");
    try appendLine(arena, &out, "");
    try appendLine(arena, &out, "[Workspace]");
    for (p.modules) |m| try appendLine(arena, &out, try arena.print("{s}=0, 0, 0, 0, C", .{m.name}));
    return out.items;
}

fn projectWmStream(arena: Allocator, p: Project) ![]u8 {
    var out: ByteArrayList = .empty;
    for (p.modules) |m| {
        try out.appendSlice(arena, m.name);
        try out.append(arena, 0);
        for (m.name) |c| try bytes.appendIntLE(arena, &out, u16, c);
        try bytes.appendIntLE(arena, &out, u16, 0);
    }
    try bytes.appendIntLE(arena, &out, u16, 0);
    return out.items;
}

/// [MS-OVBA] 2.4.3.1: the sum of the project ID's bytes, modulo 256.
fn projectKey(id: []const u8) u8 {
    var key: u8 = 0;
    for (id) |c| key +%= c;
    return key;
}

fn appendEncryptedLine(arena: Allocator, out: *ByteArrayList, comptime name: []const u8, seed: u8, key: u8, data: []const u8) !void {
    const enc = try encrypt(arena, seed, key, data);
    try out.appendSlice(arena, name ++ "=\"");
    for (enc) |b| try out.print(arena, "{X:0>2}", .{b});
    try appendLine(arena, out, "\"");
}

/// [MS-OVBA] 2.4.3.2 Data Encryption. `seed` would normally be random; a
/// fixed one keeps the output reproducible.
pub fn encrypt(gpa: Allocator, seed: u8, project_key: u8, data: []const u8) ![]u8 {
    const version: u8 = 2;
    var out: ByteArrayList = .empty;
    errdefer out.deinit(gpa);
    const version_enc = seed ^ version;
    const key_enc = seed ^ project_key;
    try out.appendSlice(gpa, &.{ seed, version_enc, key_enc });

    var unencrypted_1 = project_key;
    var encrypted_1 = key_enc;
    var encrypted_2 = version_enc;
    const State = struct {
        fn push(o: *ByteArrayList, a: Allocator, b: u8, unenc_1: *u8, enc_1: *u8, enc_2: *u8) !void {
            const enc = b ^ (enc_2.* +% unenc_1.*);
            try o.append(a, enc);
            enc_2.* = enc_1.*;
            enc_1.* = enc;
            unenc_1.* = b;
        }
    };
    const ignored_len = (seed & 6) / 2;
    for (0..ignored_len) |_| try State.push(&out, gpa, 0, &unencrypted_1, &encrypted_1, &encrypted_2);
    var len_bytes: [4]u8 = undefined;
    bytes.writeIntLE(u32, &len_bytes, @as(u32, @intCast(data.len)));
    for (len_bytes) |b| try State.push(&out, gpa, b, &unencrypted_1, &encrypted_1, &encrypted_2);
    for (data) |b| try State.push(&out, gpa, b, &unencrypted_1, &encrypted_1, &encrypted_2);
    return out.toOwnedSlice(gpa);
}

// ---------------------------------------------------------------------------
// dir stream ([MS-OVBA] 2.3.4.2), stored compressed
// ---------------------------------------------------------------------------

const DirWriter = struct {
    gpa: Allocator,
    out: ByteArrayList = .empty,

    fn int(self: *DirWriter, comptime T: type, v: T) !void {
        try bytes.appendIntLE(self.gpa, &self.out, T, v);
    }

    /// A record with a 2-byte id, 4-byte size, and `data`.
    fn record(self: *DirWriter, id: u16, data: []const u8) !void {
        try self.int(u16, id);
        try self.int(u32, @intCast(data.len));
        try self.out.appendSlice(self.gpa, data);
    }

    fn recordInt(self: *DirWriter, id: u16, comptime T: type, v: T) !void {
        try self.int(u16, id);
        try self.int(u32, @sizeOf(T));
        try self.int(T, v);
    }

    /// An ASCII string as UTF-16LE.
    fn unicode(self: *DirWriter, id: u16, s: []const u8) !void {
        try self.int(u16, id);
        try self.int(u32, @intCast(s.len * 2));
        for (s) |c| try self.int(u16, c);
    }

    /// The MBCS record followed by its Unicode twin (`reserved_id`).
    fn pair(self: *DirWriter, id: u16, reserved_id: u16, s: []const u8) !void {
        try self.record(id, s);
        try self.unicode(reserved_id, s);
    }
};

fn dirStream(arena: Allocator, p: Project) ![]u8 {
    var d: DirWriter = .{ .gpa = arena };

    // PROJECTINFORMATION
    try d.recordInt(0x0001, u32, 0x00000003); // SYSKIND: 64-bit Windows
    try d.recordInt(0x004A, u32, 0x00000006); // COMPATVERSION
    try d.recordInt(0x0002, u32, p.lcid); // LCID
    try d.recordInt(0x0014, u32, p.lcid); // LCIDINVOKE
    try d.recordInt(0x0003, u16, p.code_page); // CODEPAGE
    try d.record(0x0004, p.name); // NAME
    try d.pair(0x0005, 0x0040, ""); // DOCSTRING
    try d.pair(0x0006, 0x003D, ""); // HELPFILEPATH (path 1, then path 2)
    try d.recordInt(0x0007, u32, 0); // HELPCONTEXT
    try d.recordInt(0x0008, u32, 0); // LIBFLAGS
    // VERSION: the size field is a reserved 4, followed by major and minor.
    try d.int(u16, 0x0009);
    try d.int(u32, 4);
    try d.int(u32, 1);
    try d.int(u16, 0);
    try d.pair(0x000C, 0x003C, ""); // CONSTANTS

    // PROJECTREFERENCES
    for (p.references) |r| {
        try checkAscii(r.name);
        try checkAscii(r.libid);
        try d.pair(0x0016, 0x003E, r.name); // REFERENCENAME
        try d.int(u16, 0x000D); // REFERENCEREGISTERED
        try d.int(u32, @intCast(4 + r.libid.len + 4 + 2));
        try d.int(u32, @intCast(r.libid.len));
        try d.out.appendSlice(arena, r.libid);
        try d.int(u32, 0);
        try d.int(u16, 0);
    }

    // PROJECTMODULES
    try d.recordInt(0x000F, u16, @intCast(p.modules.len));
    try d.recordInt(0x0013, u16, 0xFFFF); // PROJECTCOOKIE
    for (p.modules) |m| {
        try d.record(0x0019, m.name); // MODULENAME
        try d.unicode(0x0047, m.name); // MODULENAMEUNICODE
        try d.pair(0x001A, 0x0032, m.name); // MODULESTREAMNAME
        try d.pair(0x001C, 0x0048, ""); // MODULEDOCSTRING
        try d.recordInt(0x0031, u32, 0); // MODULEOFFSET: no performance cache
        try d.recordInt(0x001E, u32, 0); // MODULEHELPCONTEXT
        try d.recordInt(0x002C, u16, 0xFFFF); // MODULECOOKIE
        try d.record(if (m.kind == .procedural) 0x0021 else 0x0022, ""); // MODULETYPE
        try d.record(0x002B, ""); // terminator
    }
    try d.record(0x0010, ""); // dir terminator
    return d.out.items;
}

// ---------------------------------------------------------------------------
// Compression ([MS-OVBA] 2.4.1)
// ---------------------------------------------------------------------------

const chunk_size = 4096;

const CopyTokenHelp = struct {
    /// Bits holding the offset; 4..12 for positions within a 4096-byte chunk.
    bit_count: u4,
    max_length: u16,

    /// `difference` is the current position relative to the chunk start.
    fn init(difference: usize) CopyTokenHelp {
        // Smallest integer >= log2(difference), but at least 4.
        const bits: usize = if (difference <= 1) 0 else std.math.log2_int(usize, difference - 1) + 1;
        const bit_count: u4 = @intCast(@max(bits, 4));
        return .{ .bit_count = bit_count, .max_length = (@as(u16, 0xFFFF) >> bit_count) + 3 };
    }

    fn lengthMask(self: CopyTokenHelp) u16 {
        return @as(u16, 0xFFFF) >> self.bit_count;
    }

    /// Shift that positions the offset in the high bits of a copy token.
    fn offsetShift(self: CopyTokenHelp) u4 {
        return @intCast(16 - @as(u5, self.bit_count));
    }
};

/// Compresses `src` into a CompressedContainer. Caller owns the result.
pub fn compress(gpa: Allocator, src: []const u8) ![]u8 {
    var out: ByteArrayList = .empty;
    errdefer out.deinit(gpa);
    try out.append(gpa, 0x01); // SignatureByte
    var start: usize = 0;
    while (start < src.len) {
        const end = @min(start + chunk_size, src.len);
        try compressChunk(gpa, &out, src[start..end]);
        start = end;
    }
    return out.toOwnedSlice(gpa);
}

fn compressChunk(gpa: Allocator, out: *ByteArrayList, chunk: []const u8) !void {
    const header_pos = out.items.len;
    try out.appendSlice(gpa, &.{ 0, 0 });
    var pos: usize = 0;
    while (pos < chunk.len) {
        const flag_pos = out.items.len;
        try out.append(gpa, 0);
        var flags: u8 = 0;
        for (0..8) |bit| {
            if (pos >= chunk.len) break;
            const help: CopyTokenHelp = .init(pos);
            var best_len: usize = 0;
            var best_off: usize = 0;
            var cand = pos;
            while (cand > 0) {
                cand -= 1;
                var len: usize = 0;
                while (pos + len < chunk.len and len < help.max_length and chunk[cand + len] == chunk[pos + len]) len += 1;
                if (len > best_len) {
                    best_len = len;
                    best_off = pos - cand;
                }
            }
            if (best_len >= 3) {
                const token = (@as(u16, @intCast(best_off - 1)) << help.offsetShift()) | @as(u16, @intCast(best_len - 3));
                try bytes.appendIntLE(gpa, out, u16, token);
                flags |= @as(u8, 1) << @intCast(bit);
                pos += best_len;
            } else {
                try out.append(gpa, chunk[pos]);
                pos += 1;
            }
        }
        out.items[flag_pos] = flags;
    }

    const data_len = out.items.len - header_pos - 2;
    if (data_len > chunk_size) {
        // Incompressible: store raw, which must be exactly 4096 bytes.
        if (chunk.len != chunk_size) return error.VbaChunkIncompressible;
        out.shrinkRetainingCapacity(header_pos);
        try bytes.appendIntLE(gpa, out, u16, 0x3000 | (chunk_size - 1));
        try out.appendSlice(gpa, chunk);
        return;
    }
    // Header: chunk size minus 3 (12 bits), signature 0b011, compressed flag.
    const header: u16 = 0x8000 | 0x3000 | @as(u16, @intCast(data_len + 2 - 3));
    bytes.writeIntLE(u16, out.items[header_pos..][0..2], header);
}

/// Decompresses a CompressedContainer. Caller owns the result.
pub fn decompress(gpa: Allocator, src: []const u8) ![]u8 {
    if (src.len == 0 or src[0] != 0x01) return error.VbaBadContainer;
    var out: ByteArrayList = .empty;
    errdefer out.deinit(gpa);
    var pos: usize = 1;
    while (pos < src.len) {
        if (pos + 2 > src.len) return error.VbaTruncated;
        const header = bytes.readIntLE(u16, src[pos..][0..2]);
        if ((header >> 12) & 0x7 != 0x3) return error.VbaBadChunk;
        const end = @min(pos + (header & 0x0FFF) + 3, src.len);
        pos += 2;
        const chunk_start = out.items.len;
        if (header & 0x8000 == 0) {
            if (pos + chunk_size > src.len) return error.VbaTruncated;
            try out.appendSlice(gpa, src[pos..][0..chunk_size]);
            pos += chunk_size;
            continue;
        }
        while (pos < end) {
            const flags = src[pos];
            pos += 1;
            for (0..8) |bit| {
                if (pos >= end) break;
                if ((flags >> @intCast(bit)) & 1 == 0) {
                    try out.append(gpa, src[pos]);
                    pos += 1;
                } else {
                    if (pos + 2 > end) return error.VbaTruncated;
                    const token = bytes.readIntLE(u16, src[pos..][0..2]);
                    pos += 2;
                    const help: CopyTokenHelp = .init(out.items.len - chunk_start);
                    const len = (token & help.lengthMask()) + 3;
                    const offset = (token >> help.offsetShift()) + 1;
                    if (offset > out.items.len - chunk_start) return error.VbaBadCopyToken;
                    for (0..len) |_| try out.append(gpa, out.items[out.items.len - offset]);
                }
            }
        }
    }
    return out.toOwnedSlice(gpa);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

fn expectRoundTrip(src: []const u8) !void {
    const gpa = zigstd.testAllocator;
    const c = try compress(gpa, src);
    defer gpa.free(c);
    const d = try decompress(gpa, c);
    defer gpa.free(d);
    try zigstd.testExpectEqualSlices(u8, src, d);
}

test "compression round trips" {
    try expectRoundTrip("");
    try expectRoundTrip("a");
    // cspell:disable-next-line
    try expectRoundTrip("#aaabcdefaaaaghijaaaaaklaaamnopqaaaaaaaaaaaarstuvwxyzaaa");
    var big: [10_000]u8 = undefined;
    for (&big, 0..) |*c, i| c.* = "Option Explicit\r\nPublic Sub X()\r\n"[i % 33];
    try expectRoundTrip(&big);
    var noise: [9000]u8 = undefined;
    var prng: std.Random.DefaultPrng = .init(1);
    prng.random().bytes(&noise);
    try expectRoundTrip(&noise);
}

/// Inverse of `encrypt`, for tests: returns the project key and the data.
fn decrypt(gpa: Allocator, enc: []const u8) !struct { key: u8, data: []u8 } {
    const seed = enc[0];
    if (seed ^ enc[1] != 2) return error.VbaBadVersion;
    const key = seed ^ enc[2];
    var unenc_1 = key;
    var enc_1 = enc[2];
    var enc_2 = enc[1];
    var plain: ByteArrayList = .empty;
    defer plain.deinit(gpa);
    for (enc[3..]) |e| {
        const b = e ^ (enc_2 +% unenc_1);
        try plain.append(gpa, b);
        enc_2 = enc_1;
        enc_1 = e;
        unenc_1 = b;
    }
    const ignored = (seed & 6) / 2;
    const len = bytes.readIntLE(u32, plain.items[ignored..][0..4]);
    return .{ .key = key, .data = try gpa.dupe(u8, plain.items[ignored + 4 ..][0..len]) };
}

test "encryption round trips" {
    const gpa = zigstd.testAllocator;
    for ([_]u8{ 0x00, 0x36, 0x6C, 0xA2, 0xFF }) |seed| {
        const enc = try encrypt(gpa, seed, 0xE0, &.{ 0, 0, 0, 0 });
        defer gpa.free(enc);
        const dec = try decrypt(gpa, enc);
        defer gpa.free(dec.data);
        try zigstd.testExpectEqual(@as(u8, 0xE0), dec.key);
        try zigstd.testExpectEqualSlices(u8, &.{ 0, 0, 0, 0 }, dec.data);
    }
}

test "project key is the byte sum of the ID" {
    // Excel's project {924D388A-845A-4F4F-98A7-67F259E8A661} decrypts with key 0xE0.
    try zigstd.testExpectEqual(@as(u8, 0xE0), projectKey("{924D388A-845A-4F4F-98A7-67F259E8A661}"));
}

test stripExportHeader {
    const exported =
        "VERSION 1.0 CLASS\r\nBEGIN\r\n  MultiUse = -1  'True\r\nEND\r\n" ++
        "Attribute VB_Name = \"ThisWorkbook\"\r\nAttribute VB_Exposed = True\r\n" ++
        "Option Explicit\r\n";
    try zigstd.testExpectEqualStrings("Option Explicit\r\n", stripExportHeader(exported));
    try zigstd.testExpectEqualStrings("Sub X()\n", stripExportHeader("Attribute VB_Name = \"m\"\nSub X()\n"));
}

test "module source gets attributes and CRLF" {
    var arena_state: std.heap.ArenaAllocator = .init(zigstd.testAllocator);
    defer arena_state.deinit();
    const src = try moduleSource(arena_state.allocator(), .{ .name = "shTasks", .kind = .worksheet, .code = "Option Explicit\n" });
    try zigstd.testExpectEqualStrings(
        "Attribute VB_Name = \"shTasks\"\r\n" ++
            "Attribute VB_Base = \"0{00020820-0000-0000-C000-000000000046}\"\r\n" ++
            "Attribute VB_GlobalNameSpace = False\r\n" ++
            "Attribute VB_Creatable = False\r\n" ++
            "Attribute VB_PredeclaredId = True\r\n" ++
            "Attribute VB_Exposed = True\r\n" ++
            "Attribute VB_TemplateDerived = False\r\n" ++
            "Attribute VB_Customizable = True\r\n" ++
            "Option Explicit\r\n",
        src,
    );
}
