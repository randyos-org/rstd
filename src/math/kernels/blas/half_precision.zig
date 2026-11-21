//! Which half-precision `gemv` or `gemm`: the types in and out, and the
//! family whose slots and blocking it takes.

const types = @import("types.zig");
const portable = @import("portable.zig");

pub const HalfPrecision = enum {
    /// bfloat16 in, `f32` out: `sbgemv`, `sbgemm`
    bf16_to_f32,
    /// bfloat16 throughout: `bgemv`, `bgemm`
    bf16,
    /// `f16` in, `f32` out: `shgemv`, `shgemm`
    f16_to_f32,

    /// The type of `A` and `x` or `B`.
    pub fn In(comptime k: HalfPrecision) type {
        return if (k == .f16_to_f32) types.HF16 else types.BFloat16;
    }

    /// The type of `y` or `C`, `alpha` and `beta`.
    pub fn Out(comptime k: HalfPrecision) type {
        return if (k == .bf16) types.BFloat16 else f32;
    }

    /// Its family's gemm, which names its slots and blocking parameters.
    pub fn gemm(comptime k: HalfPrecision) portable.Gemm {
        return switch (k) {
            .bf16_to_f32 => .sb,
            .bf16 => .b,
            .f16_to_f32 => .sh,
        };
    }
};
