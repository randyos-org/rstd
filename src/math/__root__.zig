//! `rstd.math`: the scalar functions here, then the kernel layer and the
//! high-level domains as they move in (`design.md`, issue #606).

const scalar = @import("scalar.zig");

pub const kernels = @import("kernels/__root__.zig");
pub const fft = @import("fft.zig");
pub const signal = @import("signal/__root__.zig");
pub const stats = @import("stats/__root__.zig");

pub const x87_f64 = scalar.x87_f64;
pub const floor = scalar.floor;
pub const ceil = scalar.ceil;
pub const round = scalar.round;
pub const Fraction = scalar.Fraction;
pub const rationalApprox = scalar.rationalApprox;

test {
    _ = scalar;
    _ = fft;
    _ = signal;
    _ = stats;
}
