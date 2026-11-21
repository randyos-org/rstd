//! [MS-OVBA] VBA project writer: builds a `vbaProject.bin` from module
//! source text, with no Office installation needed.
//!
//! Only source code is stored. `_VBA_PROJECT` carries version 0xFFFF and no
//! performance cache, which per [MS-OVBA] 2.3.4.1 tells the host there is no
//! valid compiled p-code, so Office compiles the modules from source on load.
//! Module names and source must be ASCII (a subset of every MBCS code page).
const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const bytes = @import("../../bytes.zig");
const cfb = @import("../cfb.zig");
const compression = @import("compression.zig");
const encryption = @import("encryption.zig");

const Allocator = zigstd.Allocator;
const ByteArrayList = bytes.ByteArrayList;
const projectKey = encryption.projectKey;

pub const encrypt = encryption.encrypt;
pub const compress = compression.compress;
pub const decompress = compression.decompress;

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

fn appendEncryptedLine(arena: Allocator, out: *ByteArrayList, comptime name: []const u8, seed: u8, key: u8, data: []const u8) !void {
    const enc = try encrypt(arena, seed, key, data);
    try out.appendSlice(arena, name ++ "=\"");
    for (enc) |b| try out.print(arena, "{X:0>2}", .{b});
    try appendLine(arena, out, "\"");
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
// Tests
// ---------------------------------------------------------------------------

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
