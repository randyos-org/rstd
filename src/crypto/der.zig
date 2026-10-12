//! ASN.1 DER (X.690) reading and writing, and PEM (RFC 7468) unwrapping.

const std = @import("std");
const Allocator = std.mem.Allocator;

pub const Tag = struct {
    pub const int = 0x02;
    pub const bits = 0x03;
    pub const octets = 0x04;
    pub const null_ = 0x05;
    pub const oid = 0x06;
    pub const seq = 0x30;
    pub const set = 0x31;

    /// Context-specific tag `n`, constructed (as for EXPLICIT tagging).
    pub fn context(comptime n: u5) u8 {
        return 0xA0 | @as(u8, n);
    }
};

pub const Error = error{BadDer};

pub const Tlv = struct {
    tag: u8,
    value: []const u8,
    /// Tag, length and value.
    all: []const u8,

    pub fn expect(t: Tlv, tag: u8) Error!Tlv {
        return if (t.tag == tag) t else error.BadDer;
    }

    /// INTEGER contents without the sign-padding zero.
    pub fn unsignedInt(t: Tlv) Error![]const u8 {
        const v = (try t.expect(Tag.int)).value;
        return if (v.len > 1 and v[0] == 0) v[1..] else v;
    }
};

/// Reads the element at the start of `s`.
pub fn read(s: []const u8) Error!Tlv {
    if (s.len < 2) return error.BadDer;
    var len: usize = s[1];
    var at: usize = 2;
    if (len & 0x80 != 0) {
        const n = len & 0x7F;
        if (n == 0 or n > 4 or s.len < 2 + n) return error.BadDer;
        len = 0;
        for (s[2 .. 2 + n]) |c| len = (len << 8) | c;
        at += n;
    }
    if (s.len - at < len) return error.BadDer;
    return .{ .tag = s[0], .value = s[at..][0..len], .all = s[0 .. at + len] };
}

/// Children of a constructed element's `value`, in order.
pub fn children(arena: Allocator, value: []const u8) (Error || Allocator.Error)![]Tlv {
    var out: std.ArrayList(Tlv) = .empty;
    var rest = value;
    while (rest.len > 0) {
        const t = try read(rest);
        try out.append(arena, t);
        rest = rest[t.all.len..];
    }
    return out.items;
}

/// Encodes one element whose value is the concatenation of `parts`.
pub fn encode(arena: Allocator, tag: u8, parts: []const []const u8) Allocator.Error![]u8 {
    var len: usize = 0;
    for (parts) |p| len += p.len;
    var out: std.ArrayList(u8) = .empty;
    try out.append(arena, tag);
    if (len < 0x80) {
        try out.append(arena, @intCast(len));
    } else {
        const n: u8 = if (len <= 0xFF) 1 else if (len <= 0xFFFF) 2 else if (len <= 0xFF_FFFF) 3 else 4;
        try out.append(arena, 0x80 | n);
        var i = n;
        while (i > 0) : (i -= 1) try out.append(arena, @truncate(len >> @intCast(8 * (i - 1))));
    }
    for (parts) |p| try out.appendSlice(arena, p);
    return out.items;
}

/// SET OF, with members sorted by encoding as DER requires.
pub fn encodeSetOf(arena: Allocator, members: [][]const u8) Allocator.Error![]u8 {
    std.mem.sort([]const u8, members, {}, struct {
        fn lt(_: void, a: []const u8, b: []const u8) bool {
            return std.mem.lessThan(u8, a, b);
        }
    }.lt);
    return encode(arena, Tag.set, members);
}

/// OBJECT IDENTIFIER element from its encoded arcs.
pub fn oid(comptime body: []const u8) *const [body.len + 2]u8 {
    return comptime &([2]u8{ Tag.oid, body.len } ++ body[0..body.len].*);
}

pub const null_value = [2]u8{ Tag.null_, 0 };

/// Common object identifiers, as complete elements.
pub const Oid = struct {
    pub const rsa_encryption = oid("\x2a\x86\x48\x86\xf7\x0d\x01\x01\x01");
    pub const sha256 = oid("\x60\x86\x48\x01\x65\x03\x04\x02\x01");
    pub const sha384 = oid("\x60\x86\x48\x01\x65\x03\x04\x02\x02");
    pub const sha512 = oid("\x60\x86\x48\x01\x65\x03\x04\x02\x03");
    pub const pkcs7_signed_data = oid("\x2a\x86\x48\x86\xf7\x0d\x01\x07\x02");
    pub const pkcs9_content_type = oid("\x2a\x86\x48\x86\xf7\x0d\x01\x09\x03");
    pub const pkcs9_message_digest = oid("\x2a\x86\x48\x86\xf7\x0d\x01\x09\x04");
};

/// DER body of the first PEM block labeled `label`, e.g. "CERTIFICATE".
pub fn fromPem(arena: Allocator, pem: []const u8, comptime label: []const u8) ![]u8 {
    const begin = "-----BEGIN " ++ label ++ "-----";
    const start = (std.mem.indexOf(u8, pem, begin) orelse return error.NoPemBlock) + begin.len;
    const end = std.mem.indexOfPos(u8, pem, start, "-----END ") orelse return error.BadPem;
    var b64: std.ArrayList(u8) = .empty;
    for (pem[start..end]) |c| {
        if (!std.ascii.isWhitespace(c)) try b64.append(arena, c);
    }
    const dec = std.base64.standard.Decoder;
    const out = try arena.alloc(u8, try dec.calcSizeForSlice(b64.items));
    try dec.decode(out, b64.items);
    return out;
}

test "long-form length round trip" {
    var arena_state = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena_state.deinit();
    const a = arena_state.allocator();
    const long = try a.alloc(u8, 300);
    @memset(long, 7);
    const enc = try encode(a, Tag.octets, &.{long});
    try std.testing.expectEqualSlices(u8, &.{ Tag.octets, 0x82, 0x01, 0x2C }, enc[0..4]);
    const t = try read(enc);
    try std.testing.expectEqual(@as(usize, 300), t.value.len);
}

test "set of sorts members" {
    var arena_state = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena_state.deinit();
    var members = [_][]const u8{ &.{ 0x04, 1, 9 }, &.{ 0x02, 1, 5 } };
    const enc = try encodeSetOf(arena_state.allocator(), &members);
    try std.testing.expectEqualSlices(u8, &.{ Tag.set, 6, 0x02, 1, 5, 0x04, 1, 9 }, enc);
}

test "pem" {
    var arena_state = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena_state.deinit();
    const pem = "junk\n-----BEGIN THING-----\nAQID\nBA==\n-----END THING-----\n";
    try std.testing.expectEqualSlices(u8, &.{ 1, 2, 3, 4 }, try fromPem(arena_state.allocator(), pem, "THING"));
}
