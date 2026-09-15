//! Linear algebra kernels ported from `matlab88.f90`, the 1988 Argonne
//! National Laboratory "MATLAB" Fortran source (see `matlab88_doc.txt` /
//! `matlab88_help.txt` / `readme.txt` in this directory for the original
//! documentation) - specifically its LINPACK/EISPACK-derived numerical
//! core, not its command-line interpreter.
//!
//! Each module keeps the original Fortran routine names (lowercased to
//! Zig convention; Fortran itself is case-insensitive, so this changes
//! nothing about the source) as the low-level API - `lu.wgefa`,
//! `svd.wsvdc`, `eigen_general.comqr3`, and so on. `api` layers the
//! original REPL's command names (`det`, `inv`, `eig`, `chol`, ...) on top
//! of those, including the REPL's dispatch logic (e.g. `eig` picking the
//! Hermitian vs. general eigensolver) where that logic was more than
//! cosmetic - see `api.zig`'s doc comment for what's included and why.

const std = @import("std");

pub const blas1 = @import("blas1.zig");
pub const lu = @import("lu.zig");
pub const cholesky = @import("cholesky.zig");
pub const generators = @import("generators.zig");
pub const tridiagonal = @import("tridiagonal.zig");
pub const eigen_symmetric = @import("eigen_symmetric.zig");
pub const hessenberg = @import("hessenberg.zig");
pub const eigen_general = @import("eigen_general.zig");
pub const svd = @import("svd.zig");
pub const qr = @import("qr.zig");
pub const api = @import("api.zig");

pub const Cx = blas1.Cx;

test {
    std.testing.refAllDecls(blas1);
    std.testing.refAllDecls(lu);
    std.testing.refAllDecls(cholesky);
    std.testing.refAllDecls(generators);
    std.testing.refAllDecls(tridiagonal);
    std.testing.refAllDecls(eigen_symmetric);
    std.testing.refAllDecls(hessenberg);
    std.testing.refAllDecls(eigen_general);
    std.testing.refAllDecls(svd);
    std.testing.refAllDecls(qr);
    std.testing.refAllDecls(api);
}
