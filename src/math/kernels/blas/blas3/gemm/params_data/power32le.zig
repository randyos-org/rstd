//! `power32le`'s gemm blocking parameters, from OpenBLAS `v0.3.34-794-gdbf2c24c9`'s `param.h`; converted by `tools/blas_asm` (retired: issue #618) and maintained by hand. How to read it: `__root__.zig`.

const coretypes = @import("../../../cpu/coretypes.zig");
const params = @import("../params.zig");
const Table = params.Table;
const power32 = @import("power32.zig").power32;
const power32_dtb = @import("power32.zig").power32_dtb;

/// `BUFFER_SIZE` of a core with no table of its own.
pub const power32le_buffer_size: i64 = 256 << 20;

/// `power32`'s tables but where the byte order changes them.
pub fn power32le(core: coretypes.PowerCoreType) ?Table {
    return switch (core) {
        .ppc970 => .{
            .s = .{ .unroll_m = 16, .unroll_n = 4, .p = .{ .fixed = 176 }, .q = .{ .fixed = 512 } },
            .d = .{ .unroll_m = 4, .unroll_n = 4, .p = .{ .fixed = 176 }, .q = .{ .fixed = 256 } },
            .q = .{ .unroll_m = 2, .unroll_n = 2 },
            .c = .{ .unroll_m = 8, .unroll_n = 2, .p = .{ .fixed = 176 }, .q = .{ .fixed = 256 } },
            .z = .{ .unroll_m = 2, .unroll_n = 2, .p = .{ .fixed = 176 }, .q = .{ .fixed = 128 } },
            .x = .{ .unroll_m = 2, .unroll_n = 2 },
            .sh = .{ .unroll_m = 8, .unroll_n = 8, .unroll_mn = 32, .align_k = 1, .p = .{ .fixed = 128 }, .q = .{ .fixed = 12288 }, .r = .{ .fixed = 240 } },
            .b = .{ .unroll_m = 8, .unroll_n = 4, .unroll_mn = 32, .align_k = 1, .p = .{ .fixed = 256 }, .q = .{ .fixed = 256 }, .r = .{ .fixed = 256 } },
            .sb = .{ .unroll_m = 8, .unroll_n = 4, .unroll_mn = 32, .align_k = 1, .p = .{ .fixed = 256 }, .q = .{ .fixed = 256 }, .r = .{ .fixed = 256 } },
            .offset_a = 2688,
            .offset_b = 3072,
            .align_ = 16383,
            .symv_p = 4,
            .buffer_size = 256 << 20,
        },
        else => power32(core),
    };
}

pub const power32le_dtb = power32_dtb;
