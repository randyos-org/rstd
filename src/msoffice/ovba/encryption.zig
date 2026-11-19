//! [MS-OVBA] 2.4.3 data encryption, for the `PROJECT` stream's protection fields.

const zigstd = @import("../../zigstd.zig");
const bytes = @import("../../bytes.zig");
const Allocator = zigstd.Allocator;
const ByteArrayList = bytes.ByteArrayList;

/// [MS-OVBA] 2.4.3.1: the sum of the project ID's bytes, modulo 256.
pub fn projectKey(id: []const u8) u8 {
    var key: u8 = 0;
    for (id) |c| key +%= c;
    return key;
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
