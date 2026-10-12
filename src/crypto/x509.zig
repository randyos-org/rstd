//! Just enough X.509 (RFC 5280) to sign with a certificate: issuer, serial
//! and RSA public key. std.crypto.Certificate parses for verification and
//! doesn't expose the raw issuer and serial encodings.

const std = @import("std");
const der = @import("der.zig");
const Allocator = std.mem.Allocator;

pub const Cert = struct {
    der: []const u8,
    /// Complete Name element.
    issuer: []const u8,
    /// Complete INTEGER element.
    serial: []const u8,
    /// RSA modulus of the subject key, without leading zeros.
    rsa_n: []const u8,

    pub fn fromPem(arena: Allocator, pem: []const u8) !Cert {
        return fromDer(arena, try der.fromPem(arena, pem, "CERTIFICATE"));
    }

    pub fn fromDer(arena: Allocator, bytes: []const u8) !Cert {
        const top = try der.children(arena, (try (try der.read(bytes)).expect(der.Tag.seq)).value);
        if (top.len < 1) return error.BadCert;
        var tbs = try der.children(arena, (try top[0].expect(der.Tag.seq)).value);
        if (tbs.len > 0 and tbs[0].tag == der.Tag.context(0)) tbs = tbs[1..]; // version
        // serial, signature, issuer, validity, subject, subjectPublicKeyInfo
        if (tbs.len < 6) return error.BadCert;
        const spki = try der.children(arena, (try tbs[5].expect(der.Tag.seq)).value);
        if (spki.len < 2) return error.BadCert;
        const key_bits = (try spki[1].expect(der.Tag.bits)).value;
        if (key_bits.len < 1 or key_bits[0] != 0) return error.BadCert;
        const rsa_key = try der.children(arena, (try (try der.read(key_bits[1..])).expect(der.Tag.seq)).value);
        if (rsa_key.len < 2) return error.NotRsa;
        return .{
            .der = bytes,
            .serial = (try tbs[0].expect(der.Tag.int)).all,
            .issuer = (try tbs[2].expect(der.Tag.seq)).all,
            .rsa_n = try rsa_key[0].unsignedInt(),
        };
    }
};
