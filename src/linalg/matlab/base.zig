//! Arbitrary-base number formatting, ported from `matlab88.f90`'s `BASE`.
//!
//! The original wrote its digits as indices into the REPL's own display
//! alphabet (`ALFA`/`ALFB`, used to render MATLAB "text matrices" on
//! screen) - a table this port never carried over, since it belongs to
//! the interpreter's text-matrix rendering, not to computing anything.
//! This keeps the same algorithm (normalize `|x|` to a base-`b` mantissa
//! in `[1/b, 1)` times `b^k`, with a mantissa length driven by `eps`) but
//! writes plain ASCII digits (`0`-`9`, then `A`-`Z` for bases up to 36)
//! instead, so the result is directly usable without that table.

const std = @import("std");

fn digitChar(d: usize) u8 {
    return if (d < 10) '0' + @as(u8, @intCast(d)) else 'A' + @as(u8, @intCast(d - 10));
}

/// BASE: format `x` in base `b` (`2 <= b <= 36`) as `"<sign>0.<mantissa>,<sign><exponent>"`
/// - e.g. base 10 renders like `+0.314159,+01` for `pi/10`. `eps` sets the
/// mantissa length (`-log_b(eps) + 4` digits; pass `0` for a full
/// double-precision mantissa in base `b`). Caller owns the returned slice.
pub fn base(allocator: std.mem.Allocator, x: f64, b: f64, eps: f64) ![]u8 {
    var xa = @abs(x);
    var k: i64 = 0;
    if (xa != 0.0) k = @intFromFloat(@floor(@log(xa) / @log(b)));
    if (xa > 1.0) k += 1;
    xa /= std.math.pow(f64, b, @floatFromInt(k));
    if (b * xa >= b) {
        k += 1;
        xa /= b;
    }

    const mantissa_digits: usize = if (eps != 0.0)
        @intFromFloat(@floor(-@log(eps) / @log(b) + 4.0))
    else
        54;

    var out: std.ArrayList(u8) = .empty;
    errdefer out.deinit(allocator);
    try out.append(allocator, if (x >= 0.0) '+' else '-');
    try out.appendSlice(allocator, "0.");

    var l: usize = 0;
    while (l < mantissa_digits) : (l += 1) {
        xa = b * xa;
        const j: usize = @intFromFloat(@floor(xa));
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
            const j: usize = @intFromFloat(@mod(t, b));
            try exp_digits.append(allocator, digitChar(j));
            t = @floor(t / b);
        }
    }
    std.mem.reverse(u8, exp_digits.items);
    try out.appendSlice(allocator, exp_digits.items);
    return out.toOwnedSlice(allocator);
}

const testing = std.testing;

test "base renders a simple base-10 value" {
    // 3.14 = 0.314 * 10^1.
    const s = try base(testing.allocator, 3.14, 10.0, 0.001);
    defer testing.allocator.free(s);
    try testing.expect(std.mem.startsWith(u8, s, "+0.314"));
    try testing.expect(std.mem.endsWith(u8, s, ",+1"));
}

test "base marks negative numbers and zero exponent" {
    // -0.5 = -0.5 * 10^0.
    const s = try base(testing.allocator, -0.5, 10.0, 0.01);
    defer testing.allocator.free(s);
    try testing.expect(std.mem.startsWith(u8, s, "-0.5"));
    try testing.expect(std.mem.endsWith(u8, s, ",+0"));
}
