//! A `CFF2` table (OpenType 1.9.1, "CFF2"), read into a `CffTable`: its
//! header and Top DICT, the CFF2 INDEXes with 32-bit counts, the FDArray
//! every CFF2 font has, an optional FDSelect, the variation store's offset
//! and `maxstack` (FreeType 2.14.3 `src/cff/cffload.c` `cff_font_load`).

const std = @import("std");
const bytes = @import("../../../../../bytes.zig");
const CffTable = @import("CffTable.zig");
const Index = @import("Index.zig");
const Dict = @import("Dict.zig");

/// FreeType's `CFF2_DEFAULT_STACK`.
pub const default_maxstack = 513;

pub fn init(table_data: []const u8, check: bool) (bytes.BytesError || CffTable.CffError)!CffTable {
    if (table_data.len < 5 or table_data[0] != 2) return error.EndOfStream;
    const header_size: usize = table_data[2];
    const top_len: usize = std.mem.readInt(u16, table_data[3..5], .big);
    if (header_size + top_len > table_data.len) return error.EndOfStream;
    const top: Dict = .{ .table_data = table_data[header_size..][0..top_len] };
    const global_subrs_index = try Index.initCff2(table_data[header_size + top_len ..], check);

    const charstrings_index = try Index.initCff2(try at(table_data, try top.getInt(.char_strings) orelse return error.MissingCharStrings), check);
    const fd_array_index = try Index.initCff2(try at(table_data, try top.getInt(.fd_array) orelse return error.MissingFdArray), check);
    const fdselect_offset: ?usize = if (try top.getInt(.fd_select)) |o| try offset(table_data, o) else null;
    const vstore_offset: ?usize = if (try top.getInt(.vstore)) |o| try offset(table_data, o) else null;
    const maxstack = try top.getInt(.maxstack) orelse default_maxstack;
    return .{
        .table_data = table_data,
        .top_dict_data = top.table_data,
        .charstrings_index = charstrings_index,
        .global_subrs_index = global_subrs_index,
        .is_cid = false,
        .local_subrs_index = null,
        .fd_array_index = fd_array_index,
        .fdselect_offset = fdselect_offset,
        .cff2 = true,
        .vstore_offset = vstore_offset,
        .maxstack = @intCast(std.math.clamp(maxstack, 0, 0xFFFF)),
    };
}

fn offset(table_data: []const u8, raw: i32) bytes.BytesError!usize {
    if (raw < 0 or raw >= table_data.len) return error.EndOfStream;
    return @intCast(raw);
}

fn at(table_data: []const u8, raw: i32) bytes.BytesError![]const u8 {
    return table_data[try offset(table_data, raw)..];
}
