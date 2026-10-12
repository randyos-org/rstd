//! RSA private keys and RSASSA-PKCS1-v1_5 signing (RFC 8017). std has RSA
//! verification only.

const std = @import("std");
const der = @import("der.zig");
const Allocator = std.mem.Allocator;

pub const PrivateKey = struct {
    /// Big-endian, without leading zeros.
    n: []const u8,
    e: []const u8,
    d: []const u8,

    /// PKCS#8 ("PRIVATE KEY", as openssl writes) or PKCS#1 ("RSA PRIVATE KEY") PEM.
    pub fn fromPem(arena: Allocator, pem: []const u8) !PrivateKey {
        if (der.fromPem(arena, pem, "PRIVATE KEY")) |pkcs8| {
            const info = try der.children(arena, (try (try der.read(pkcs8)).expect(der.Tag.seq)).value);
            if (info.len < 3) return error.BadKey;
            return fromPkcs1(arena, (try info[2].expect(der.Tag.octets)).value);
        } else |_| return fromPkcs1(arena, try der.fromPem(arena, pem, "RSA PRIVATE KEY"));
    }

    /// RSAPrivateKey DER.
    pub fn fromPkcs1(arena: Allocator, pkcs1: []const u8) !PrivateKey {
        const k = try der.children(arena, (try (try der.read(pkcs1)).expect(der.Tag.seq)).value);
        if (k.len < 4) return error.BadKey;
        return .{ .n = try k[1].unsignedInt(), .e = try k[2].unsignedInt(), .d = try k[3].unsignedInt() };
    }
};

/// Signature of a `Hash` digest; `Hash` is SHA-256, -384 or -512.
pub fn signPkcs1v15(arena: Allocator, key: PrivateKey, comptime Hash: type, digest: [Hash.digest_length]u8) ![]u8 {
    const M = std.crypto.ff.Modulus(4096);
    const hash_oid = switch (Hash) {
        std.crypto.hash.sha2.Sha256 => der.Oid.sha256,
        std.crypto.hash.sha2.Sha384 => der.Oid.sha384,
        std.crypto.hash.sha2.Sha512 => der.Oid.sha512,
        else => @compileError("unsupported hash"),
    };
    // DigestInfo ::= SEQUENCE { AlgorithmIdentifier, OCTET STRING digest }
    const alg_len = hash_oid.len + der.null_value.len;
    const prefix = [_]u8{ der.Tag.seq, alg_len + 4 + digest.len, der.Tag.seq, alg_len } ++
        hash_oid.* ++ der.null_value ++ [_]u8{ der.Tag.octets, digest.len };

    const k = key.n.len;
    if (k < prefix.len + digest.len + 11) return error.KeyTooSmall;
    // EM = 00 01 FF..FF 00 DigestInfo
    const em = try arena.alloc(u8, k);
    const t_at = k - prefix.len - digest.len;
    em[0] = 0;
    em[1] = 1;
    @memset(em[2 .. t_at - 1], 0xFF);
    em[t_at - 1] = 0;
    @memcpy(em[t_at..][0..prefix.len], &prefix);
    @memcpy(em[k - digest.len ..], &digest);

    const m = try M.fromBytes(key.n, .big);
    const sig = try m.powWithEncodedExponent(try M.Fe.fromBytes(m, em, .big), key.d, .big);
    const out = try arena.alloc(u8, k);
    try sig.toBytes(out, .big);

    // Catches a d that doesn't belong with n and e.
    const back = try arena.alloc(u8, k);
    try (try m.powWithEncodedPublicExponent(sig, key.e, .big)).toBytes(back, .big);
    if (!std.mem.eql(u8, back, em)) return error.InconsistentKey;
    return out;
}

test "sign with a textbook key" {
    var arena_state = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena_state.deinit();
    // Too small for a DigestInfo; exercises the size check only.
    const key: PrivateKey = .{ .n = &.{ 0x0C, 0xA1 }, .e = &.{0x11}, .d = &.{ 0x0A, 0xC1 } };
    var digest: [32]u8 = undefined;
    @memset(&digest, 0);
    try std.testing.expectError(error.KeyTooSmall, signPkcs1v15(arena_state.allocator(), key, std.crypto.hash.sha2.Sha256, digest));
}
