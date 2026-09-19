//! Compile-time helpers for the transcribed `.S` kernels, which are asm
//! templates written as the upstream source with its macros left as names.

const std = @import("std");

/// One macro: every whole-word `name` in a template becomes `text`.
pub const Sub = struct { name: []const u8, text: []const u8 };

/// `template` with the `subs` applied, `%` doubled for Zig's asm syntax, and
/// each local `.L` label suffixed with `suffix`, so kernels that reuse label
/// names can share an object.
pub fn expand(comptime template: []const u8, comptime subs: []const Sub, comptime suffix: []const u8) []const u8 {
    comptime {
        @setEvalBranchQuota(1 << 26);
        const text = substitute(template, subs);
        var buf: [text.len * 2 + text.len / 4 * suffix.len]u8 = undefined;
        var n: usize = 0;
        var i: usize = 0;
        while (i < text.len) {
            if (text[i] == '%') {
                buf[n..][0..2].* = "%%".*;
                n += 2;
                i += 1;
            } else if (text[i] == '.' and i + 1 < text.len and text[i + 1] == 'L' and (i == 0 or !isIdent(text[i - 1]))) {
                var j = i + 2;
                while (j < text.len and isIdent(text[j])) j += 1;
                @memcpy(buf[n..][0 .. j - i], text[i..j]);
                n += j - i;
                @memcpy(buf[n..][0..suffix.len], suffix);
                n += suffix.len;
                i = j;
            } else {
                buf[n] = text[i];
                n += 1;
                i += 1;
            }
        }
        const out: [n]u8 = buf[0..n].*;
        return &out;
    }
}

fn substitute(comptime template: []const u8, comptime subs: []const Sub) []const u8 {
    comptime {
        var longest: usize = 1;
        for (subs) |s| longest = @max(longest, s.text.len);
        var buf: [template.len * longest]u8 = undefined;
        var n: usize = 0;
        var i: usize = 0;
        while (i < template.len) {
            const c = template[i];
            // numbers (`0x3e`, `1f`) and words are taken whole
            if (isIdent(c)) {
                var j = i;
                while (j < template.len and isIdent(template[j])) j += 1;
                const word = template[i..j];
                const text = if (std.ascii.isDigit(c)) word else for (subs) |s| {
                    if (std.mem.eql(u8, s.name, word)) break s.text;
                } else word;
                @memcpy(buf[n..][0..text.len], text);
                n += text.len;
                i = j;
            } else {
                buf[n] = c;
                n += 1;
                i += 1;
            }
        }
        const out: [n]u8 = buf[0..n].*;
        return &out;
    }
}

fn isIdent(c: u8) bool {
    return std.ascii.isAlphanumeric(c) or c == '_';
}
