//! Prints the CPU and the BLAS core rstd detects and runs on it; CI's test
//! legs show it before their tests (issue #619).

const std = @import("std");
const sysinfo = @import("builtin");
const rstd = @import("rstd");

const blas = rstd.linalg.blas;

pub fn main() void {
    std.debug.print("{t}-{t}, cpu model {s}, detected core {t}, runs {s}\n", .{
        sysinfo.target.cpu.arch,
        sysinfo.target.os.tag,
        sysinfo.target.cpu.model.name,
        blas.cpu.detect(),
        blas.coreName(),
    });
}
