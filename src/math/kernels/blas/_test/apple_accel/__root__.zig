//! `apple_accel.zig`'s tests: binding, routing, the split over rstd's
//! workers, and special values against rstd's kernels.

test {
    _ = @import("binding.zig");
    _ = @import("routing.zig");
    _ = @import("splits.zig");
    _ = @import("special.zig");
}
