//! The x86_64 dot drivers (`blas1/dot/driver.zig`) against OpenBLAS's own C
//! (`kernel/x86_64/sdot.c`, `ddot.c`, `cdot.c`, `zdot.c`), bit for bit, fused
//! and not, on the cores this CPU runs. Not a test: `tools/dotverify/run.sh`
//! builds it with the C, under a stub `common.h`, and runs it.

const std = @import("std");
const sysinfo = @import("builtin");
const kernel = @import("blas1/dot/driver.zig");
const x86 = @import("blas1/dot/x86_64/kernel.zig");
const ukidx = @import("blas1/dot/x86_64/uk/__root__.zig");
const util = @import("util.zig");

const Slot = util.Slot;
const Feature = std.Target.x86.Feature;

/// The cores `run.sh` builds the C for, with what running them needs here.
const cores = [_]struct { name: []const u8, needs: []const Feature }{
    .{ .name = "nehalem", .needs = &.{.sse4_1} },
    .{ .name = "sandybridge", .needs = &.{.avx} },
    .{ .name = "haswell", .needs = &.{ .avx2, .fma } },
    .{ .name = "zen", .needs = &.{ .avx2, .fma } },
    .{ .name = "skylakex", .needs = &.{.avx2} },
    .{ .name = "generic", .needs = &.{} },
};

const slots = [_]Slot{ .dsdot_k, .ddot_k, .cdotu_k, .cdotc_k, .zdotu_k, .zdotc_k };

const sizes = [_]usize{ 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 34, 63, 64, 65, 100, 127, 128, 129, 257, 640, 1279, 1280, 1281, 1500, 4100 };
const strides = [_][2]usize{ .{ 1, 1 }, .{ 2, 1 }, .{ 1, 3 }, .{ 2, 2 } };

/// The reference's block size where a core has no block of its own.
fn fallbackMask(comptime slot: Slot) comptime_int {
    return switch (slot) {
        .dsdot_k => 32,
        .ddot_k, .cdotu_k, .cdotc_k => 16,
        else => 8,
    };
}

/// `slot`'s driver on `core`, as `uk_block.driven` builds it, with `fuse`.
fn Driver(comptime slot: Slot, comptime core: []const u8, comptime fuse: bool) type {
    const e = ukidx.find(slot, core) orelse return kernel.Driver(slot, fallbackMask(slot), kernel.Fallback(slot, fuse), fuse);
    const B = switch (e.block) {
        .none => kernel.Fallback(slot, fuse),
        .vector => x86.VectorBlockWith(slot, std.Target.x86.featureSetHas(sysinfo.target.cpu.features, .avx512cd), fuse),
        .assembly => kernel.Asm(slot, @field(e.module, e.kernel)),
    };
    return kernel.Driver(slot, e.mask, B, fuse);
}

fn fill(comptime T: type, buf: []T, seed: u64) void {
    var prng = std.Random.DefaultPrng.init(seed);
    const r = prng.random();
    for (buf) |*v| v.* = @floatCast(r.float(f64) * 4 - 2);
}

var runs: usize = 0;
var fails: usize = 0;

fn check(comptime slot: Slot, comptime core: []const u8, comptime fuse: bool) void {
    const E = kernel.form(slot).E;
    const sym = std.fmt.comptimePrint("ref_{s}_{s}_{s}", .{ @tagName(slot)[0 .. @tagName(slot).len - 2], core, if (fuse) "f" else "n" });
    const ref = @extern(kernel.Kernel(slot), .{ .name = sym });
    const D = Driver(slot, core, fuse);
    const S = struct {
        var x: [sizes[sizes.len - 1] * 3 * 2]E = undefined;
        var y: [sizes[sizes.len - 1] * 3 * 2]E = undefined;
    };
    for (sizes) |n| {
        fill(E, &S.x, n + 1);
        fill(E, &S.y, n + 7);
        for (strides) |inc| {
            const want = ref(@intCast(n), &S.x, @intCast(inc[0]), &S.y, @intCast(inc[1]));
            const got = D.call(@intCast(n), &S.x, @intCast(inc[0]), &S.y, @intCast(inc[1]));
            runs += 1;
            if (!std.mem.eql(u8, std.mem.asBytes(&want), std.mem.asBytes(&got))) {
                fails += 1;
                if (fails <= 40) std.debug.print("MISMATCH {t} {s} fuse={} n={d} inc={d},{d}: ref {any} got {any}\n", .{ slot, core, fuse, n, inc[0], inc[1], want, got });
            }
        }
    }
}

pub fn main() u8 {
    inline for (cores) |c| {
        const runnable = comptime for (c.needs) |f| {
            if (!std.Target.x86.featureSetHas(sysinfo.target.cpu.features, f)) break false;
        } else true;
        if (runnable) {
            inline for (.{ false, true }) |fuse| inline for (slots) |slot| check(slot, c.name, fuse);
        } else std.debug.print("{s}: not run, this CPU lacks its features\n", .{c.name});
    }
    std.debug.print("{d} runs, {d} mismatches\n", .{ runs, fails });
    return @intFromBool(fails != 0);
}
