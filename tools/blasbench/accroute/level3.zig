//! The `accroute` section's level 3: `syr2k`, `hemm`, `herk`, `her2k`, square
//! n, upper, no transpose, `alpha` and `beta` 1 (`hemm` from the left).

const rstd = @import("rstd");
const bb = @import("../blas_bench.zig");
const libs = @import("../libs.zig");
const aw = @import("work.zig");

const blas = rstd.math.kernels.blas;
const Prec = bb.Prec;
const rstd_alloc = libs.rstd_alloc;
const Work = aw.Work;
const scalar = aw.scalar;

pub const cases = .{ Syr2k, Hemm, Herk, Her2k };

const Syr2k = struct {
    pub const name = "syr2k";
    pub const letters = "sdcz";
    pub const shape: aw.Shape = .mm;
    pub const fill: aw.Fill = .general;
    pub const restore: aw.Restore = .none;
    pub const out: aw.Out = .y;
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        const one = scalar(p, 1, 0);
        try blas.syr2k(rstd_alloc, libs.rstd_io, Work(p).T, .u, .n, w.n, w.n, one, w.a, w.n, w.b, w.n, one, w.y, w.n);
    }
};

const Hemm = struct {
    pub const name = "hemm";
    pub const letters = "cz";
    pub const shape: aw.Shape = .mm;
    pub const fill: aw.Fill = .general;
    pub const restore: aw.Restore = .none;
    pub const out: aw.Out = .y;
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        const one = scalar(p, 1, 0);
        try blas.hemm(rstd_alloc, libs.rstd_io, Work(p).T, .l, .u, w.n, w.n, one, w.a, w.n, w.b, w.n, one, w.y, w.n);
    }
};

const Herk = struct {
    pub const name = "herk";
    pub const letters = "cz";
    pub const shape: aw.Shape = .mm;
    pub const fill: aw.Fill = .general;
    pub const restore: aw.Restore = .none;
    pub const out: aw.Out = .y;
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.herk(rstd_alloc, libs.rstd_io, Work(p).T, .u, .n, w.n, w.n, 1, w.a, w.n, 1, w.y, w.n);
    }
};

const Her2k = struct {
    pub const name = "her2k";
    pub const letters = "cz";
    pub const shape: aw.Shape = .mm;
    pub const fill: aw.Fill = .general;
    pub const restore: aw.Restore = .none;
    pub const out: aw.Out = .y;
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.her2k(rstd_alloc, libs.rstd_io, Work(p).T, .u, .n, w.n, w.n, scalar(p, 1, 0), w.a, w.n, w.b, w.n, 1, w.y, w.n);
    }
};
