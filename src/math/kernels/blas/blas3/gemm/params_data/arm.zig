//! `arm`'s gemm blocking parameters, from OpenBLAS `v0.3.34-794-gdbf2c24c9`'s `param.h`; converted by `tools/blas_asm` (retired: issue #618) and maintained by hand. How to read it: `__root__.zig`.

const coretypes = @import("../../../cpu/coretypes.zig");
const params = @import("../params.zig");
const Table = params.Table;

/// `BUFFER_SIZE` of a core with no table of its own.
pub const arm_buffer_size: i64 = 32 << 20;

pub fn arm(core: coretypes.ArmCoreType) ?Table {
    return switch (core) {
        .armv5 => .{
            .s = .{ .unroll_m = 2, .unroll_n = 2, .p = .{ .fixed = 128 }, .q = .{ .fixed = 240 }, .r = .{ .fixed = 12288 } },
            .d = .{ .unroll_m = 2, .unroll_n = 2, .p = .{ .fixed = 128 }, .q = .{ .fixed = 120 }, .r = .{ .fixed = 8192 } },
            .q = .{ .unroll_m = 2, .unroll_n = 2 },
            .c = .{ .unroll_m = 2, .unroll_n = 2, .p = .{ .fixed = 96 }, .q = .{ .fixed = 120 }, .r = .{ .fixed = 4096 } },
            .z = .{ .unroll_m = 2, .unroll_n = 2, .p = .{ .fixed = 64 }, .q = .{ .fixed = 120 }, .r = .{ .fixed = 4096 } },
            .x = .{ .unroll_m = 2, .unroll_n = 2 },
            .sh = .{ .unroll_m = 8, .unroll_n = 8, .unroll_mn = 32, .align_k = 1, .p = .{ .fixed = 128 }, .q = .{ .fixed = 12288 }, .r = .{ .fixed = 240 } },
            .b = .{ .unroll_m = 8, .unroll_n = 4, .unroll_mn = 32, .align_k = 1, .p = .{ .fixed = 256 }, .q = .{ .fixed = 256 }, .r = .{ .fixed = 256 } },
            .sb = .{ .unroll_m = 8, .unroll_n = 4, .unroll_mn = 32, .align_k = 1, .p = .{ .fixed = 256 }, .q = .{ .fixed = 256 }, .r = .{ .fixed = 256 } },
            .align_ = 16383,
            .symv_p = 16,
            .buffer_size = 32 << 20,
        },
        .armv6 => .{
            .s = .{ .unroll_m = 4, .unroll_n = 2, .p = .{ .fixed = 128 }, .q = .{ .fixed = 240 }, .r = .{ .fixed = 12288 } },
            .d = .{ .unroll_m = 4, .unroll_n = 2, .p = .{ .fixed = 128 }, .q = .{ .fixed = 120 }, .r = .{ .fixed = 8192 } },
            .q = .{ .unroll_m = 2, .unroll_n = 2 },
            .c = .{ .unroll_m = 2, .unroll_n = 2, .p = .{ .fixed = 96 }, .q = .{ .fixed = 120 }, .r = .{ .fixed = 4096 } },
            .z = .{ .unroll_m = 2, .unroll_n = 2, .p = .{ .fixed = 64 }, .q = .{ .fixed = 120 }, .r = .{ .fixed = 4096 } },
            .x = .{ .unroll_m = 2, .unroll_n = 2 },
            .sh = .{ .unroll_m = 8, .unroll_n = 8, .unroll_mn = 32, .align_k = 1, .p = .{ .fixed = 128 }, .q = .{ .fixed = 12288 }, .r = .{ .fixed = 240 } },
            .b = .{ .unroll_m = 8, .unroll_n = 4, .unroll_mn = 32, .align_k = 1, .p = .{ .fixed = 256 }, .q = .{ .fixed = 256 }, .r = .{ .fixed = 256 } },
            .sb = .{ .unroll_m = 8, .unroll_n = 4, .unroll_mn = 32, .align_k = 1, .p = .{ .fixed = 256 }, .q = .{ .fixed = 256 }, .r = .{ .fixed = 256 } },
            .align_ = 16383,
            .symv_p = 16,
            .buffer_size = 32 << 20,
        },
        .armv7, .neon, .vfpv4, .cortexa9, .cortexa15 => .{
            .s = .{ .unroll_m = 4, .unroll_n = 4, .p = .{ .fixed = 128 }, .q = .{ .fixed = 240 }, .r = .{ .fixed = 12288 } },
            .d = .{ .unroll_m = 4, .unroll_n = 4, .p = .{ .fixed = 128 }, .q = .{ .fixed = 120 }, .r = .{ .fixed = 8192 } },
            .q = .{ .unroll_m = 2, .unroll_n = 2 },
            .c = .{ .unroll_m = 2, .unroll_n = 2, .p = .{ .fixed = 96 }, .q = .{ .fixed = 120 }, .r = .{ .fixed = 4096 } },
            .z = .{ .unroll_m = 2, .unroll_n = 2, .p = .{ .fixed = 64 }, .q = .{ .fixed = 120 }, .r = .{ .fixed = 4096 } },
            .x = .{ .unroll_m = 2, .unroll_n = 2 },
            .sh = .{ .unroll_m = 8, .unroll_n = 8, .unroll_mn = 32, .align_k = 1, .p = .{ .fixed = 128 }, .q = .{ .fixed = 12288 }, .r = .{ .fixed = 240 } },
            .b = .{ .unroll_m = 8, .unroll_n = 4, .unroll_mn = 32, .align_k = 1, .p = .{ .fixed = 256 }, .q = .{ .fixed = 256 }, .r = .{ .fixed = 256 } },
            .sb = .{ .unroll_m = 8, .unroll_n = 4, .unroll_mn = 32, .align_k = 1, .p = .{ .fixed = 256 }, .q = .{ .fixed = 256 }, .r = .{ .fixed = 256 } },
            .align_ = 16383,
            .symv_p = 16,
            .buffer_size = 32 << 20,
        },
    };
}

/// `DTB_DEFAULT_ENTRIES`; null where `getarch.c` has no entry for the core.
pub fn arm_dtb(core: coretypes.ArmCoreType) ?i32 {
    return switch (core) {
        .armv5, .armv6, .armv7, .neon => 64,
        .vfpv4, .cortexa9, .cortexa15 => 128,
    };
}
