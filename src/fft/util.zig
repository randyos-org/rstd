//! PocketFFT's size helpers (`pocketfft_hdronly.h`'s `util`, cpp branch,
//! commit 4d865cba, BSD-3, LICENSE): the plan's choice between mixed radix
//! and Bluestein's.

const std = @import("std");

pub fn largestPrimeFactor(n0: usize) usize {
    var n = n0;
    var res: usize = 1;
    while (n & 1 == 0) {
        res = 2;
        n >>= 1;
    }
    var x: usize = 3;
    while (x * x <= n) : (x += 2) {
        while (n % x == 0) {
            res = x;
            n /= x;
        }
    }
    if (n > 1) res = n;
    return res;
}

/// The cost model: factors up to 5 cost themselves, larger ones 1.1 times.
pub fn costGuess(n0: usize) f64 {
    const lfp = 1.1; // penalty for non-hardcoded larger factors
    var n = n0;
    var result: f64 = 0;
    while (n & 1 == 0) {
        result += 2;
        n >>= 1;
    }
    var x: usize = 3;
    while (x * x <= n) : (x += 2) {
        while (n % x == 0) {
            const fx: f64 = @floatFromInt(x);
            result += if (x <= 5) fx else lfp * fx;
            n /= x;
        }
    }
    if (n > 1) {
        const fnn: f64 = @floatFromInt(n);
        result += if (n <= 5) fnn else lfp * fnn;
    }
    return result * @as(f64, @floatFromInt(n0));
}

/// The smallest product of 2, 3, 5, 7 and 11 that is at least `n`; null
/// past PocketFFT's overflow guard.
pub fn goodSizeComplex(n: usize) ?usize {
    if (n <= 12) return n;
    if (n > std.math.maxInt(usize) / 11 / 2) return null;
    var bestfac = 2 * n;
    var f11: usize = 1;
    while (f11 < bestfac) : (f11 *= 11) {
        var f117 = f11;
        while (f117 < bestfac) : (f117 *= 7) {
            var f1175 = f117;
            while (f1175 < bestfac) : (f1175 *= 5) {
                var x = f1175;
                while (x < n) x *= 2;
                while (true) {
                    if (x < n) {
                        x *= 3;
                    } else if (x > n) {
                        if (x < bestfac) bestfac = x;
                        if (x & 1 != 0) break;
                        x >>= 1;
                    } else return n;
                }
            }
        }
    }
    return bestfac;
}

test "goodSizeComplex: the next 2, 3, 5, 7, 11-smooth size" {
    const testing = std.testing;
    try testing.expectEqual(@as(?usize, 1), goodSizeComplex(1));
    try testing.expectEqual(@as(?usize, 14), goodSizeComplex(13));
    try testing.expectEqual(@as(?usize, 2025), goodSizeComplex(2017));
    try testing.expectEqual(@as(?usize, 1024), goodSizeComplex(1024));
    try testing.expectEqual(@as(usize, 1009), largestPrimeFactor(1009));
    try testing.expectEqual(@as(usize, 5), largestPrimeFactor(360));
}
