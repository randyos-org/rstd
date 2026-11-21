//! `rstd.math.kernels.stats` (`../../design.md`): summary-statistics
//! reductions, MKL's VS; the deviation sums so far (#605).

const deviation = @import("deviation.zig");

pub const DeviationSums = deviation.DeviationSums;
pub const deviationSums = deviation.deviationSums;

test {
    _ = @import("_test/deviation.zig");
}
