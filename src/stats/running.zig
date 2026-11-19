//! Streaming mean and variance.

const moments = @import("moments.zig");
const util = @import("util.zig");

const Normalize = moments.Normalize;
const checkFloat = util.checkFloat;

/// Streaming mean and variance (Welford's update). Its bits differ from
/// `mean`'s and `variance`'s.
pub fn Running(comptime F: type) type {
    comptime checkFloat(F);
    return struct {
        count: u64 = 0,
        avg: F = 0,
        /// Sum of squared deviations from `avg`.
        m2: F = 0,

        const Self = @This();

        pub fn push(self: *Self, x: F) void {
            self.count += 1;
            const d = x - self.avg;
            self.avg += d / @as(F, @floatFromInt(self.count));
            self.m2 += d * (x - self.avg);
        }

        pub fn mean(self: Self) error{EmptyInput}!F {
            if (self.count == 0) return error.EmptyInput;
            return self.avg;
        }

        pub fn variance(self: Self, normalize: Normalize) error{EmptyInput}!F {
            if (self.count == 0) return error.EmptyInput;
            const d = switch (normalize) {
                .sample => if (self.count == 1) return 0 else self.count - 1,
                .population => self.count,
            };
            return self.m2 / @as(F, @floatFromInt(d));
        }
    };
}
