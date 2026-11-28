//! The `accroute` section's level 1: `rotm`, `dsdot`, `sdsdot`, `axpby`.

const rstd = @import("rstd");
const bb = @import("../blas_bench.zig");
const libs = @import("../libs.zig");
const aw = @import("work.zig");

const blas = rstd.math.kernels.blas;
const Prec = bb.Prec;
const Work = aw.Work;
const scalar = aw.scalar;

pub const cases = .{ Rotm, Dsdot, Sdsdot, Axpby };

/// A rotation by 0.3 radians as `rotm`'s full `H` (flag -1), so the
/// vectors keep their size.
const Rotm = struct {
    pub const name = "rotm";
    pub const letters = "sd";
    pub const shape: aw.Shape = .vec;
    pub const fill: aw.Fill = .general;
    pub const restore: aw.Restore = .none;
    pub const out: aw.Out = .y;
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        const c: p.F = @cos(@as(p.F, 0.3));
        const s: p.F = @sin(@as(p.F, 0.3));
        const param = [5]p.F{ -1, c, -s, s, c };
        try blas.rotm(libs.rstd_io, p.F, w.xv(), w.yv(), &param);
    }
};

const Dsdot = struct {
    pub const name = "dsdot";
    pub const letters = "s";
    pub const shape: aw.Shape = .vec;
    pub const fill: aw.Fill = .general;
    pub const restore: aw.Restore = .none;
    pub const out: aw.Out = .r;
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        w.r = @floatCast(try blas.dsdot(libs.rstd_io, w.xc(), w.yc()));
    }
};

const Sdsdot = struct {
    pub const name = "sdsdot";
    pub const letters = "s";
    pub const shape: aw.Shape = .vec;
    pub const fill: aw.Fill = .general;
    pub const restore: aw.Restore = .none;
    pub const out: aw.Out = .r;
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        w.r = try blas.sdsdot(libs.rstd_io, 0.25, w.xc(), w.yc());
    }
};

/// `y := alpha x + beta y` with both scalars about one half: `y` stays bounded.
const Axpby = struct {
    pub const name = "axpby";
    pub const letters = "sdcz";
    pub const shape: aw.Shape = .vec;
    pub const fill: aw.Fill = .general;
    pub const restore: aw.Restore = .none;
    pub const out: aw.Out = .y;
    pub fn call(comptime p: Prec, w: *Work(p)) !void {
        try blas.axpby(libs.rstd_io, Work(p).T, scalar(p, 0.5, 0.25), w.xc(), scalar(p, 0.5, -0.25), w.yv());
    }
};
