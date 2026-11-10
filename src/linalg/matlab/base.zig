//! BASE, from `matlab88.f90`'s `BASE`, writing ASCII digits in place of the
//! REPL's display alphabet.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const math = @import("../../math.zig");
const Allocator = std.mem.Allocator;

fn digitChar(d: usize) u8 {
    return if (d < 10) '0' + @as(u8, @intCast(d)) else 'A' + @as(u8, @intCast(d - 10));
}

/// BASE: `x` in base `b` (2 to 36) as `"<sign>0.<mantissa>,<sign><exponent>"`;
/// `-log_b(eps) + 4` digits, or 54 for `eps` 0. Caller frees.
pub fn base(allocator: Allocator, x: f64, b: f64, eps: f64) ![]u8 {
    var xa = @abs(x);
    var k: i64 = 0;
    if (xa != 0.0) k = @trunc(math.floor(@log(xa) / @log(b)));
    if (xa > 1.0) k += 1;
    xa /= std.math.pow(f64, b, @floatFromInt(k));
    if (b * xa >= b) {
        k += 1;
        xa /= b;
    }

    const mantissa_digits: usize = if (eps != 0.0)
        @trunc(math.floor(-@log(eps) / @log(b) + 4.0))
    else
        54;

    var out: std.ArrayList(u8) = .empty;
    errdefer out.deinit(allocator);
    try out.append(allocator, if (x >= 0.0) '+' else '-');
    try out.appendSlice(allocator, "0.");

    var l: usize = 0;
    while (l < mantissa_digits) : (l += 1) {
        xa = b * xa;
        const j: usize = @trunc(math.floor(xa));
        try out.append(allocator, digitChar(j));
        xa -= @floatFromInt(j);
    }
    try out.append(allocator, ',');
    try out.append(allocator, if (k >= 0) '+' else '-');

    var t = @abs(@as(f64, @floatFromInt(k)));
    var exp_digits: std.ArrayList(u8) = .empty;
    defer exp_digits.deinit(allocator);
    if (t == 0.0) {
        try exp_digits.append(allocator, '0');
    } else {
        while (t > 0.0) {
            const j: usize = @trunc(@mod(t, b));
            try exp_digits.append(allocator, digitChar(j));
            t = math.floor(t / b);
        }
    }
    std.mem.reverse(u8, exp_digits.items);
    try out.appendSlice(allocator, exp_digits.items);
    return out.toOwnedSlice(allocator);
}

// ---- tests ------------------------------------------------------------------

test "base renders a simple base-10 value" {
    // 3.14 = 0.314 * 10^1.
    const s = try base(zigstd.testAllocator, 3.14, 10.0, 0.001);
    defer zigstd.testAllocator.free(s);
    try zigstd.testExpectTrue(std.mem.startsWith(u8, s, "+0.314"));
    try zigstd.testExpectTrue(std.mem.endsWith(u8, s, ",+1"));
}

test "base marks negative numbers and zero exponent" {
    // -0.5 = -0.5 * 10^0.
    const s = try base(zigstd.testAllocator, -0.5, 10.0, 0.01);
    defer zigstd.testAllocator.free(s);
    try zigstd.testExpectTrue(std.mem.startsWith(u8, s, "-0.5"));
    try zigstd.testExpectTrue(std.mem.endsWith(u8, s, ",+0"));
}
