//! Descriptive statistics over slices (rstd's design.md, "Statistics"): the
//! moments through the BLAS, `Running`, the median and quantiles,
//! histograms, and `cov` and `corrcoef`.
//!
//! NaN propagates, as in MATLAB. `cov` and `corrcoef` take the public
//! layer's views; the forms along a dimension come with N-d arrays.

const moments = @import("moments.zig");
const running = @import("running.zig");
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

pub const Running = running.Running;

pub const MedianError = quantiles.MedianError;
pub const median = quantiles.median;
pub const QuantileMethod = quantiles.QuantileMethod;
pub const QuantileOptions = quantiles.QuantileOptions;
pub const QuantileError = quantiles.QuantileError;
pub const QuantileIntoError = quantiles.QuantileIntoError;
pub const quantile = quantiles.quantile;
pub const quantileInto = quantiles.quantileInto;

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
