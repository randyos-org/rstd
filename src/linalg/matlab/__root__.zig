//! Zig port of `matlab88.f90`'s function library - the 1988 Argonne
//! National Laboratory "MATLAB" Fortran source (see `ref/matlab88_doc.txt`
//! / `ref/matlab88_help.txt` / `ref/readme.txt` for the original
//! documentation) - everything callable as a MATLAB function, minus the
//! REPL itself (parser, tokenizer, `<>` stack machine, `FOR`/`WHILE`/`IF`,
//! `SAVE`/`LOAD`/`DIARY`/`PLOT` I/O, `HELP` text) since those have no
//! meaning outside an interactive session. See `README.md` for a usage
//! guide and a full function table.
//!
//! Each module keeps the original Fortran routine names (lowercased to
//! Zig convention; Fortran itself is case-insensitive, so this changes
//! nothing about the source) as the low-level API - `lu.wgefa`,
//! `svd.wsvdc`, `eigen_general.comqr3`, and so on. `api` layers the
//! original REPL's command names (`det`, `inv`, `eig`, `chol`, `kron`,
//! `roots`, ...) on top of those, including the REPL's dispatch logic
//! (e.g. `eig` picking the Hermitian vs. general eigensolver) where that
//! logic was more than cosmetic - see `api.zig`'s doc comment for what's
//! included and why.

const std = @import("std");

pub const blas1 = @import("blas1.zig");
pub const elementary = @import("elementary.zig");
pub const lu = @import("lu.zig");
pub const cholesky = @import("cholesky.zig");
pub const generators = @import("generators.zig");
pub const tridiagonal = @import("tridiagonal.zig");
pub const eigen_symmetric = @import("eigen_symmetric.zig");
pub const hessenberg = @import("hessenberg.zig");
pub const eigen_general = @import("eigen_general.zig");
pub const svd = @import("svd.zig");
pub const qr = @import("qr.zig");
pub const rat = @import("rat.zig");
pub const base = @import("base.zig");
pub const random = @import("random.zig");
pub const eig = @import("eig.zig");
pub const divide = @import("divide.zig");
pub const svdapi = @import("svdapi.zig");
pub const matrixfun = @import("matrixfun.zig");
pub const matmul = @import("matmul.zig");
pub const schur = @import("schur.zig");
pub const norm = @import("norm.zig");
pub const api = @import("api.zig");

pub const Cx = blas1.Cx;

test {
    std.testing.refAllDecls(blas1);
    std.testing.refAllDecls(elementary);
    std.testing.refAllDecls(lu);
    std.testing.refAllDecls(cholesky);
    std.testing.refAllDecls(generators);
    std.testing.refAllDecls(tridiagonal);
    std.testing.refAllDecls(eigen_symmetric);
    std.testing.refAllDecls(hessenberg);
    std.testing.refAllDecls(eigen_general);
    std.testing.refAllDecls(svd);
    std.testing.refAllDecls(qr);
    std.testing.refAllDecls(rat);
    std.testing.refAllDecls(base);
    std.testing.refAllDecls(random);
    std.testing.refAllDecls(eig);
    std.testing.refAllDecls(divide);
    std.testing.refAllDecls(svdapi);
    std.testing.refAllDecls(matrixfun);
    std.testing.refAllDecls(matmul);
    std.testing.refAllDecls(schur);
    std.testing.refAllDecls(norm);
    std.testing.refAllDecls(api);
}
