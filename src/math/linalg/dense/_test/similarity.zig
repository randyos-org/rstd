//! Test helper: the residual of a similarity transform.

const matrix = @import("../../matrix.zig");

const Mat = matrix.Mat;

/// `‖Q X Qᵀ - A‖` over the entries, real `A`.
pub fn residual(a: []const f64, q: Mat(f64), x: Mat(f64)) f64 {
    const n = q.rows;
    var worst: f64 = 0;
    for (0..n) |i| for (0..n) |j| {
        var s: f64 = 0;
        for (0..n) |k| for (0..n) |l| {
            s += q.at(i, k) * x.at(k, l) * q.at(j, l);
        };
        worst = @max(worst, @abs(s - a[i + j * n]));
    };
    return worst;
}
