const std = @import("std");
const uefi = std.os.uefi;
const Io = std.Io;
const log = std.log.scoped(.random);

const state = @import("state.zig");
const time = @import("time.zig");

/// Outer optional = "never looked" vs "looked, absent"; locateProtocol is a
/// firmware call, only do it once. Module global not `Threaded` field: only
/// ever one RNG protocol.
var rng_lookup: ??*uefi.protocol.Rng = null;

fn locateRng() ?*uefi.protocol.Rng {
    if (rng_lookup == null) {
        const bs = state.bootServices() orelse return null;
        rng_lookup = bs.locateProtocol(uefi.protocol.Rng, null) catch null;
    }
    return rng_lookup.?;
}

// ------------------------------
// io vtable callbacks below here

pub fn randomSecure(userdata: ?*anyopaque, buffer: []u8) Io.RandomSecureError!void {
    _ = userdata;
    const rng = locateRng() orelse return error.EntropyUnavailable;
    rng.getRNG(null, buffer) catch return error.EntropyUnavailable;
}

/// Fallback for `random` when firmware has no RNG protocol; seeded from tick
/// clock on first use. Not crypto-secure -- use `randomSecure` for that.
var fallback_prng: ?std.Random.DefaultPrng = null;

pub fn random(userdata: ?*anyopaque, buffer: []u8) void {
    randomSecure(userdata, buffer) catch {
        if (fallback_prng == null) {
            var seed: u64 = time.getElapsedNanoseconds();
            // stack addr adds per-boot variety when timer hasn't started (seed=0)
            seed ^= @intFromPtr(&seed) *% 0x9e3779b97f4a7c15;
            fallback_prng = .init(seed);
        }
        fallback_prng.?.fill(buffer);
    };
}
