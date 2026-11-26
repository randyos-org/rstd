//! `rstd.math.kernels.stats` (`../../design.md`): summary-statistics
//! kernels, MKL's VS: the deviation sums (#605), the running moments and
//! histogram counts.

const deviation = @import("deviation.zig");
const running = @import("running.zig");
const histogram = @import("histogram.zig");

pub const DeviationSums = deviation.DeviationSums;
pub const deviationSums = deviation.deviationSums;

pub const Normalize = running.Normalize;
pub const RunningMoments = running.RunningMoments;

pub const Bins = histogram.Bins;
pub const Histogram = histogram.Histogram;
pub const HistogramIntoError = histogram.HistogramIntoError;
pub const histogramInto = histogram.histogramInto;
pub const checkBins = histogram.checkBins;
pub const bin = histogram.bin;

test {
    _ = @import("_test/deviation.zig");
    _ = @import("_test/histogram.zig");
}
