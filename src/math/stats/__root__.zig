//! `rstd.math.stats` (`../design.md`, #613): descriptive statistics over
//! slices: the moments through the kernel layer, `RunningMoments`, the
//! median and quantiles, histograms, and `cov` and `corrcoef`.
//!
//! NaN propagates, as in MATLAB. `cov` and `corrcoef` take the public
//! layer's views; the forms along a dimension come with N-d arrays.

const kernels_stats = @import("../kernels/stats/__root__.zig");
const moments = @import("moments.zig");
const quantiles = @import("quantile.zig");
const histograms = @import("histogram.zig");
const covariance = @import("cov.zig");

pub const Normalize = moments.Normalize;
pub const MomentError = moments.MomentError;
pub const Summation = moments.Summation;
pub const MomentOptions = moments.MomentOptions;
pub const mean = moments.mean;
pub const meanWith = moments.meanWith;
pub const variance = moments.variance;
pub const varianceWith = moments.varianceWith;
pub const stdDev = moments.stdDev;
pub const stdDevWith = moments.stdDevWith;

pub const RunningMoments = kernels_stats.RunningMoments;

pub const MedianError = quantiles.MedianError;
pub const median = quantiles.median;
pub const QuantileMethod = quantiles.QuantileMethod;
pub const QuantileOptions = quantiles.QuantileOptions;
pub const QuantileError = quantiles.QuantileError;
pub const QuantileIntoError = quantiles.QuantileIntoError;
pub const quantile = quantiles.quantile;
pub const quantileWith = quantiles.quantileWith;
pub const quantileInto = quantiles.quantileInto;
pub const quantileIntoWith = quantiles.quantileIntoWith;

pub const Bins = histograms.Bins;
pub const Histogram = histograms.Histogram;
pub const HistogramError = histograms.HistogramError;
pub const HistogramIntoError = histograms.HistogramIntoError;
pub const histogram = histograms.histogram;
pub const histogramInto = histograms.histogramInto;

pub const CovError = covariance.CovError;
pub const cov = covariance.cov;
pub const corrcoef = covariance.corrcoef;

test {
    _ = @import("_test/stats.zig");
}
