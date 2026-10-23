//! RISC-V `axpy` on the vector extension: OpenBLAS `kernel/riscv64/axpy_rvv.c`
//! and `zaxpy_rvv.c`, for the cores whose `KERNEL` file names them
//! (`RISCV64_ZVL128B`, `x280`).
//!
//! The C is a `vsetvl` loop around vector intrinsics. Zig has no intrinsics
//! for the vector extension, so each loop is the assembly the intrinsics
//! stand for, with the same registers grouping (`m8` for real, `m4` for the
//! complex segments), the same branches on the strides, and the same
//! `vl`-driven pointer stepping.
//!
//! The loop of each branch is one `asm` block, since the vector length the
//! `vsetvli` sets has to survive to the loads and stores that use it.

const std = @import("std");
const util = @import("../../../util.zig");

/// Element width in bits and its log2 in bytes.
fn sew(comptime T: type) []const u8 {
    return if (T == f64) "64" else "32";
}

fn shift(comptime T: type) []const u8 {
    return if (T == f64) "3" else "2";
}

/// The vector registers the blocks use, and what they alias.
const clobbers: std.lang.assembly.Clobbers = .{
    .memory = true,
    .vl = true,
    .vtype = true,
    .t0 = true,
    .t1 = true,
    .t2 = true,
    .ft0 = true,
    .v0 = true,
    .v8 = true,
    .v9 = true,
    .v10 = true,
    .v11 = true,
    .v12 = true,
    .v13 = true,
    .v14 = true,
    .v15 = true,
    .v16 = true,
    .v17 = true,
    .v18 = true,
    .v19 = true,
    .v20 = true,
    .v21 = true,
    .v22 = true,
    .v23 = true,
};

/// `y := alpha*x + y`, contiguous or strided, the reference's kernel ABI.
pub fn Real(comptime T: type) type {
    const e = sew(T);
    const s = shift(T);
    const bytes = @sizeOf(T);
    return struct {
        pub fn call(
            n: isize,
            _: isize,
            _: isize,
            alpha: T,
            x: [*]T,
            inc_x: isize,
            y: [*]T,
            inc_y: isize,
            _: [*]T,
            _: isize,
        ) callconv(.c) i32 {
            if (n <= 0) return 0;
            var m: usize = @intCast(n);
            var xp = x;
            var yp = y;
            const sx: isize = inc_x * bytes;
            const sy: isize = inc_y * bytes;

            if (inc_x == 1 and inc_y == 1) {
                asm volatile (std.fmt.comptimePrint(
                        \\1:
                        \\ vsetvli t0, %[n], e{0s}, m8, ta, ma
                        \\ vle{0s}.v v8, (%[x])
                        \\ vle{0s}.v v16, (%[y])
                        \\ vfmacc.vf v16, %[a], v8
                        \\ vse{0s}.v v16, (%[y])
                        \\ sub %[n], %[n], t0
                        \\ slli t0, t0, {1s}
                        \\ add %[x], %[x], t0
                        \\ add %[y], %[y], t0
                        \\ bnez %[n], 1b
                    , .{ e, s })
                    : [n] "+r" (m),
                      [x] "+r" (xp),
                      [y] "+r" (yp),
                    : [a] "f" (alpha),
                    : clobbers);
            } else if (inc_y == 1) {
                asm volatile (std.fmt.comptimePrint(
                        \\1:
                        \\ vsetvli t0, %[n], e{0s}, m8, ta, ma
                        \\ vlse{0s}.v v8, (%[x]), %[sx]
                        \\ vle{0s}.v v16, (%[y])
                        \\ vfmacc.vf v16, %[a], v8
                        \\ vse{0s}.v v16, (%[y])
                        \\ sub %[n], %[n], t0
                        \\ mul t1, t0, %[sx]
                        \\ add %[x], %[x], t1
                        \\ slli t0, t0, {1s}
                        \\ add %[y], %[y], t0
                        \\ bnez %[n], 1b
                    , .{ e, s })
                    : [n] "+r" (m),
                      [x] "+r" (xp),
                      [y] "+r" (yp),
                    : [a] "f" (alpha),
                      [sx] "r" (sx),
                    : clobbers);
            } else if (inc_x == 1 and inc_y != 0) {
                asm volatile (std.fmt.comptimePrint(
                        \\1:
                        \\ vsetvli t0, %[n], e{0s}, m8, ta, ma
                        \\ vle{0s}.v v8, (%[x])
                        \\ vlse{0s}.v v16, (%[y]), %[sy]
                        \\ vfmacc.vf v16, %[a], v8
                        \\ vsse{0s}.v v16, (%[y]), %[sy]
                        \\ sub %[n], %[n], t0
                        \\ mul t1, t0, %[sy]
                        \\ add %[y], %[y], t1
                        \\ slli t0, t0, {1s}
                        \\ add %[x], %[x], t0
                        \\ bnez %[n], 1b
                    , .{ e, s })
                    : [n] "+r" (m),
                      [x] "+r" (xp),
                      [y] "+r" (yp),
                    : [a] "f" (alpha),
                      [sy] "r" (sy),
                    : clobbers);
            } else if (inc_y == 0) {
                // one at a time: the reference started every lane at y[0]
                // and summed them, so it added y[0] once a lane
                var i: isize = 0;
                while (i < n) : (i += 1) y[0] = @mulAdd(T, alpha, util.at(T, x, i * inc_x).*, y[0]);
            } else {
                asm volatile (std.fmt.comptimePrint(
                        \\1:
                        \\ vsetvli t0, %[n], e{0s}, m8, ta, ma
                        \\ vlse{0s}.v v8, (%[x]), %[sx]
                        \\ vlse{0s}.v v16, (%[y]), %[sy]
                        \\ vfmacc.vf v16, %[a], v8
                        \\ vsse{0s}.v v16, (%[y]), %[sy]
                        \\ sub %[n], %[n], t0
                        \\ mul t1, t0, %[sx]
                        \\ add %[x], %[x], t1
                        \\ mul t1, t0, %[sy]
                        \\ add %[y], %[y], t1
                        \\ bnez %[n], 1b
                    , .{e})
                    : [n] "+r" (m),
                      [x] "+r" (xp),
                      [y] "+r" (yp),
                    : [a] "f" (alpha),
                      [sx] "r" (sx),
                      [sy] "r" (sy),
                    : clobbers);
            }
            return 0;
        }
    };
}

/// The complex kernel; `conj` selects the `axpyc` formula.
pub fn Cplx(comptime F: type, comptime conj: bool) type {
    const e = sew(F);
    // a complex element is two floats: the shift of its byte size
    const s = if (F == f64) "4" else "3";
    const bytes = 2 * @sizeOf(F);
    // vy0 += ar*vx0 [- ai*vx1]; vy1 += [-]ar*vx1 + ai*vx0, as the C orders them
    const body = if (!conj)
        \\ vfmacc.vf v16, %[ar], v8
        \\ vfnmsac.vf v16, %[ai], v12
        \\ vfmacc.vf v20, %[ar], v12
        \\ vfmacc.vf v20, %[ai], v8
    else
        \\ vfmacc.vf v16, %[ar], v8
        \\ vfmacc.vf v16, %[ai], v12
        \\ vfnmsac.vf v20, %[ar], v12
        \\ vfmacc.vf v20, %[ai], v8
    ;
    return struct {
        pub fn call(
            n: isize,
            _: isize,
            _: isize,
            alpha_r: F,
            alpha_i: F,
            x: [*]F,
            inc_x: isize,
            y: [*]F,
            inc_y: isize,
            _: [*]F,
            _: isize,
        ) callconv(.c) i32 {
            if (n < 0) return 0;
            var m: usize = @intCast(n);
            var xp = x;
            var yp = y;
            const sx: isize = inc_x * bytes;
            const sy: isize = inc_y * bytes;
            if (m == 0) return 0;

            if (inc_y == 0) {
                // one at a time, each lane's operations in turn: a segment
                // loads `y` before storing it
                var i: isize = 0;
                while (i < n) : (i += 1) {
                    const xr = util.at(F, x, 2 * i * inc_x).*;
                    const xi = util.at(F, x, 2 * i * inc_x + 1).*;
                    if (!conj) {
                        y[0] = @mulAdd(F, -alpha_i, xi, @mulAdd(F, alpha_r, xr, y[0]));
                        y[1] = @mulAdd(F, alpha_i, xr, @mulAdd(F, alpha_r, xi, y[1]));
                    } else {
                        y[0] = @mulAdd(F, alpha_i, xi, @mulAdd(F, alpha_r, xr, y[0]));
                        y[1] = @mulAdd(F, alpha_i, xr, @mulAdd(F, -alpha_r, xi, y[1]));
                    }
                }
            } else if (inc_x == 1 and inc_y == 1) {
                asm volatile (std.fmt.comptimePrint(
                        \\1:
                        \\ vsetvli t0, %[n], e{0s}, m4, ta, ma
                        \\ vlseg2e{0s}.v v8, (%[x])
                        \\ vlseg2e{0s}.v v16, (%[y])
                        \\{2s}
                        \\ vsseg2e{0s}.v v16, (%[y])
                        \\ sub %[n], %[n], t0
                        \\ slli t0, t0, {1s}
                        \\ add %[x], %[x], t0
                        \\ add %[y], %[y], t0
                        \\ bnez %[n], 1b
                    , .{ e, s, body })
                    : [n] "+r" (m),
                      [x] "+r" (xp),
                      [y] "+r" (yp),
                    : [ar] "f" (alpha_r),
                      [ai] "f" (alpha_i),
                    : clobbers);
            } else if (inc_x == 1) {
                asm volatile (std.fmt.comptimePrint(
                        \\1:
                        \\ vsetvli t0, %[n], e{0s}, m4, ta, ma
                        \\ vlseg2e{0s}.v v8, (%[x])
                        \\ vlsseg2e{0s}.v v16, (%[y]), %[sy]
                        \\{2s}
                        \\ vssseg2e{0s}.v v16, (%[y]), %[sy]
                        \\ sub %[n], %[n], t0
                        \\ mul t1, t0, %[sy]
                        \\ add %[y], %[y], t1
                        \\ slli t0, t0, {1s}
                        \\ add %[x], %[x], t0
                        \\ bnez %[n], 1b
                    , .{ e, s, body })
                    : [n] "+r" (m),
                      [x] "+r" (xp),
                      [y] "+r" (yp),
                    : [ar] "f" (alpha_r),
                      [ai] "f" (alpha_i),
                      [sy] "r" (sy),
                    : clobbers);
            } else if (inc_y == 1) {
                asm volatile (std.fmt.comptimePrint(
                        \\1:
                        \\ vsetvli t0, %[n], e{0s}, m4, ta, ma
                        \\ vlsseg2e{0s}.v v8, (%[x]), %[sx]
                        \\ vlseg2e{0s}.v v16, (%[y])
                        \\{2s}
                        \\ vsseg2e{0s}.v v16, (%[y])
                        \\ sub %[n], %[n], t0
                        \\ mul t1, t0, %[sx]
                        \\ add %[x], %[x], t1
                        \\ slli t0, t0, {1s}
                        \\ add %[y], %[y], t0
                        \\ bnez %[n], 1b
                    , .{ e, s, body })
                    : [n] "+r" (m),
                      [x] "+r" (xp),
                      [y] "+r" (yp),
                    : [ar] "f" (alpha_r),
                      [ai] "f" (alpha_i),
                      [sx] "r" (sx),
                    : clobbers);
            } else {
                asm volatile (std.fmt.comptimePrint(
                        \\1:
                        \\ vsetvli t0, %[n], e{0s}, m4, ta, ma
                        \\ vlsseg2e{0s}.v v8, (%[x]), %[sx]
                        \\ vlsseg2e{0s}.v v16, (%[y]), %[sy]
                        \\{1s}
                        \\ vssseg2e{0s}.v v16, (%[y]), %[sy]
                        \\ sub %[n], %[n], t0
                        \\ mul t1, t0, %[sx]
                        \\ add %[x], %[x], t1
                        \\ mul t1, t0, %[sy]
                        \\ add %[y], %[y], t1
                        \\ bnez %[n], 1b
                    , .{ e, body })
                    : [n] "+r" (m),
                      [x] "+r" (xp),
                      [y] "+r" (yp),
                    : [ar] "f" (alpha_r),
                      [ai] "f" (alpha_i),
                      [sx] "r" (sx),
                      [sy] "r" (sy),
                    : clobbers);
            }
            return 0;
        }
    };
}
