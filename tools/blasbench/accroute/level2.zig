//! The `accroute` section's level 2: the symmetric and Hermitian products and
//! rank updates (full, band and packed) and the triangular products and solves
//! the sweep of issue #772 did not time; upper, no transpose, non-unit. The
//! products add `1e-3 A x` to `y`, the updates `1e-6` of theirs to `A`.

const rstd = @import("rstd");
const bb = @import("../blas_bench.zig");
const libs = @import("../libs.zig");
const aw = @import("work.zig");

const blas = rstd.math.kernels.blas;
const Prec = bb.Prec;
const rstd_alloc = libs.rstd_alloc;
const Work = aw.Work;
const band_k = aw.band_k;
const band_rows = aw.band_rows;
const scalar = aw.scalar;

pub const cases = .{ Hemv, Hbmv, Hpmv, Sbmv, Spmv, Her, Hpr, Her2, Hpr2, Syr, Spr, Syr2, Spr2, Trmv, Tbmv, Tpmv, Tbsv, Tpsv };

/// A case: `C.call(p, w)` makes the call; a triangular one (`out` `x`)
/// starts each call from a fresh `x`.
fn Case(comptime name_: []const u8, comptime letters_: []const u8, comptime shape_: aw.Shape, comptime fill_: aw.Fill, comptime out_: aw.Out, comptime C: type) type {
    return struct {
        pub const name = name_;
        pub const letters = letters_;
        pub const shape = shape_;
        pub const fill = fill_;
        pub const out = out_;
        pub const restore: aw.Restore = if (out_ == .x) .x else .none;
        pub const call = C.call;
    };
}

fn T(comptime p: Prec) type {
    return Work(p).T;
}
fn mvAlpha(comptime p: Prec) T(p) {
    return scalar(p, 1e-3, 5e-4);
}
fn rankAlpha(comptime p: Prec) T(p) {
    return scalar(p, 1e-6, 5e-7);
}
fn one(comptime p: Prec) T(p) {
    return scalar(p, 1, 0);
}

const Hemv = Case("hemv", "cz", .full, .general, .y, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.hemv(rstd_alloc, libs.rstd_io, T(p), .u, mvAlpha(p), w.a, w.n, w.xc(), one(p), w.yv());
    }
});
const Hbmv = Case("hbmv", "cz", .band, .general, .y, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.hbmv(rstd_alloc, libs.rstd_io, T(p), .u, band_k, mvAlpha(p), w.a, band_rows, w.xc(), one(p), w.yv());
    }
});
const Hpmv = Case("hpmv", "cz", .packed_, .general, .y, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.hpmv(rstd_alloc, libs.rstd_io, T(p), .u, mvAlpha(p), w.a, w.xc(), one(p), w.yv());
    }
});
const Sbmv = Case("sbmv", "sd", .band, .general, .y, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.sbmv(rstd_alloc, libs.rstd_io, T(p), .u, band_k, mvAlpha(p), w.a, band_rows, w.xc(), one(p), w.yv());
    }
});
const Spmv = Case("spmv", "sd", .packed_, .general, .y, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.spmv(rstd_alloc, libs.rstd_io, T(p), .u, mvAlpha(p), w.a, w.xc(), one(p), w.yv());
    }
});
const Her = Case("her", "cz", .full, .general, .a, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.her(rstd_alloc, libs.rstd_io, T(p), .u, 1e-6, w.xc(), w.a, w.n);
    }
});
const Hpr = Case("hpr", "cz", .packed_, .general, .a, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.hpr(rstd_alloc, libs.rstd_io, T(p), .u, 1e-6, w.xc(), w.a);
    }
});
const Her2 = Case("her2", "cz", .full, .general, .a, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.her2(rstd_alloc, libs.rstd_io, T(p), .u, rankAlpha(p), w.xc(), w.yc(), w.a, w.n);
    }
});
const Hpr2 = Case("hpr2", "cz", .packed_, .general, .a, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.hpr2(rstd_alloc, libs.rstd_io, T(p), .u, rankAlpha(p), w.xc(), w.yc(), w.a);
    }
});
const Syr = Case("syr", "sd", .full, .general, .a, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.syr(rstd_alloc, libs.rstd_io, T(p), .u, rankAlpha(p), w.xc(), w.a, w.n);
    }
});
const Spr = Case("spr", "sd", .packed_, .general, .a, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.spr(rstd_alloc, libs.rstd_io, T(p), .u, rankAlpha(p), w.xc(), w.a);
    }
});
const Syr2 = Case("syr2", "sd", .full, .general, .a, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.syr2(rstd_alloc, libs.rstd_io, T(p), .u, rankAlpha(p), w.xc(), w.yc(), w.a, w.n);
    }
});
const Spr2 = Case("spr2", "sd", .packed_, .general, .a, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.spr2(rstd_alloc, libs.rstd_io, T(p), .u, rankAlpha(p), w.xc(), w.yc(), w.a);
    }
});
const Trmv = Case("trmv", "sdcz", .full, .general, .x, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.trmv(rstd_alloc, libs.rstd_io, T(p), .u, .n, .n, w.a, w.n, w.xv());
    }
});
const Tbmv = Case("tbmv", "sdcz", .band, .general, .x, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.tbmv(rstd_alloc, libs.rstd_io, T(p), .u, .n, .n, band_k, w.a, band_rows, w.xv());
    }
});
const Tpmv = Case("tpmv", "sdcz", .packed_, .general, .x, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.tpmv(rstd_alloc, libs.rstd_io, T(p), .u, .n, .n, w.a, w.xv());
    }
});
const Tbsv = Case("tbsv", "sdcz", .band, .dominant, .x, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.tbsv(rstd_alloc, libs.rstd_io, T(p), .u, .n, .n, band_k, w.a, band_rows, w.xv());
    }
});
const Tpsv = Case("tpsv", "sdcz", .packed_, .dominant, .x, struct {
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.tpsv(rstd_alloc, libs.rstd_io, T(p), .u, .n, .n, w.a, w.xv());
    }
});
