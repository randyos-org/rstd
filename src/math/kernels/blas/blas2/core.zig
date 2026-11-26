//! Level 2's reference-shaped core in one place (pointers, sizes, strides,
//! the reference's named errors), for the LAPACK routines, the ported tests and
//! the other internal callers; `api.zig` is the public layer over it.

const iface = @import("interface.zig");
const tri = @import("triangular.zig");
const band = @import("banded.zig");
const rank = @import("rank.zig");

pub const Error = iface.Error;
pub const HalfError = iface.HalfError;
pub const ArgError = iface.ArgError;
pub const HalfPrecision = iface.HalfPrecision;

pub const gemv = iface.gemv;
pub const gemvHalf = iface.gemvHalf;
pub const ger = iface.ger;
pub const gerc = iface.gerc;
pub const symv = iface.symv;
pub const hemv = iface.hemv;
pub const trmv = tri.trmv;
pub const trsv = tri.trsv;
pub const tbmv = tri.tbmv;
pub const tbsv = tri.tbsv;
pub const tpmv = tri.tpmv;
pub const tpsv = tri.tpsv;
pub const gbmv = band.gbmv;
pub const sbmv = band.sbmv;
pub const hbmv = band.hbmv;
pub const spmv = band.spmv;
pub const hpmv = band.hpmv;
pub const syr = rank.syr;
pub const spr = rank.spr;
pub const her = rank.her;
pub const hpr = rank.hpr;
pub const syr2 = rank.syr2;
pub const spr2 = rank.spr2;
pub const her2 = rank.her2;
pub const hpr2 = rank.hpr2;
