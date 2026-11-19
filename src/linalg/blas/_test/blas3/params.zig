//! `blas3/gemm/params.zig`'s tests.

const std = @import("std");
const Vtable = @import("../../Vtable.zig");
const data = params.tables;
const coretypes = @import("../../cpu/coretypes.zig");
const corespec = @import("../../cpu/corespec.zig");
const select = @import("../../cpu/select.zig");
const driver = @import("../../blas3/driver.zig");
const slots = @import("../../slots.zig");
const params = @import("../../blas3/gemm/params.zig");

const Gemm = slots.Gemm;
const Param = slots.Param;
const Env = params.Env;
const detect = params.detect;
const derivedR = params.derivedR;
const Family = params.Family;
const CoreOf = params.CoreOf;
const tableFor = params.tableFor;
const dtbFor = params.dtbFor;
const GenericBranch = params.GenericBranch;
const genericBranch = params.genericBranch;
const installGeneric = params.installGeneric;
const genericBufferSize = params.genericBufferSize;
const native_buffer_size = params.native_buffer_size;
const installCore = params.installCore;

test "dtb entries: getarch's per core, cpuid's or generic's otherwise" {
    try std.testing.expectEqual(@as(i32, 256), dtbFor(.x86_64, .core2));
    try std.testing.expectEqual(@as(i32, 40), dtbFor(.x86_64, .bobcat));
    // `FORCE_P2`'s, whose `CORENAME` is `P5`
    try std.testing.expectEqual(@as(i32, 64), dtbFor(.x86, .p5));
    try std.testing.expectEqual(@as(i32, 64), dtbFor(.riscv64, .x280));
    try std.testing.expectEqual(@as(i32, 128), dtbFor(.x86_64, .unknown));
    try std.testing.expectEqual(@as(i32, 128), dtbFor(.generic, .generic));
}

/// One core's checks, in a frame of its own: a vtable is large, and the test
/// visits hundreds of cores.
noinline fn checkCore(comptime family: Family, comptime core: CoreOf(family)) !void {
    var vt: Vtable = .{};
    installCore(&vt, family, core, .{ .l2_kb = 512, .l3_kb = 8192, .threads = 1 });
    if (vt.symv_p != 4 and vt.symv_p != 8 and vt.symv_p != 16) {
        std.debug.print("{s}/{s}: symv_p={d}\n", .{ @tagName(family), @tagName(core), vt.symv_p });
        return error.TestUnexpectedResult;
    }
    inline for ([_]Gemm{ .s, .d, .c, .z }) |g| {
        const p = slots.paramPtr(&vt, g, .p).*;
        const q = slots.paramPtr(&vt, g, .q).*;
        const r = slots.paramPtr(&vt, g, .r).*;
        const um = slots.paramPtr(&vt, g, .unroll_m).*;
        // a precision the core has no header set for stays zero; only x86
        // rounds `p` to the unroll
        const x86 = family == .x86 or family == .x86_64;
        if (um != 0 and (p <= 0 or q <= 0 or r <= 0 or (x86 and @rem(p, um) != 0))) {
            std.debug.print("{s}/{s}/{s}: p={d} q={d} r={d} unroll_m={d}\n", .{ @tagName(family), @tagName(core), @tagName(g), p, q, r, um });
            return error.TestUnexpectedResult;
        }
    }
    // `Q` follows the core and caches, never the thread count (LA464's sets
    // may change `P` and `R`; threading.md, Guardrails)
    for ([_]i32{ 8192, 32 * 1024 }) |l3| {
        installCore(&vt, family, core, .{ .l2_kb = 512, .l3_kb = l3, .threads = 1 });
        const qs = qOf(&vt);
        for ([_]i32{ 0, 2, 3, 64 }) |threads| {
            installCore(&vt, family, core, .{ .l2_kb = 512, .l3_kb = l3, .threads = threads });
            if (!std.mem.eql(i32, &qs, &qOf(&vt))) {
                std.debug.print("{s}/{s}: q {any} at 1 thread, {any} at {d}\n", .{ @tagName(family), @tagName(core), qs, qOf(&vt), threads });
                return error.TestUnexpectedResult;
            }
        }
    }
}

const q_gemms = [_]Gemm{ .s, .d, .c, .z, .q, .x, .c3m, .z3m, .x3m, .sh, .b, .sb };

/// Every precision's `q`.
fn qOf(vt: *const Vtable) [q_gemms.len]i32 {
    var out: [q_gemms.len]i32 = undefined;
    inline for (q_gemms, 0..) |g, i| out[i] = slots.paramPtr(vt, g, .q).*;
    return out;
}

/// Every core of `family` resolves to sane parameters; one test per family,
/// so each fits a rig image (`build/qemu/machines.zig`'s `sparc64-debug`).
fn checkFamily(comptime family: Family) !void {
    @setEvalBranchQuota(2_000_000);
    inline for (@typeInfo(CoreOf(family)).@"enum".field_names) |cname| {
        try checkCore(family, @field(CoreOf(family), cname));
    }
}

test "every core of the x86_64 family resolves to sane parameters" {
    try checkFamily(.x86_64);
}

test "every core of the x86 family resolves to sane parameters" {
    try checkFamily(.x86);
}

test "every core of the arm64 family resolves to sane parameters" {
    try checkFamily(.arm64);
}

test "every core of the arm family resolves to sane parameters" {
    try checkFamily(.arm);
}

test "every core of the power64 family resolves to sane parameters" {
    try checkFamily(.power64);
}

test "every core of the power64be family resolves to sane parameters" {
    try checkFamily(.power64be);
}

test "every core of the power32 family resolves to sane parameters" {
    try checkFamily(.power32);
}

test "every core of the power32le family resolves to sane parameters" {
    try checkFamily(.power32le);
}

test "every core of the riscv64 family resolves to sane parameters" {
    try checkFamily(.riscv64);
}

test "every core of the s390x family resolves to sane parameters" {
    try checkFamily(.s390x);
}

test "every core of the mips64 family resolves to sane parameters" {
    try checkFamily(.mips64);
}

test "every core of the mips32 family resolves to sane parameters" {
    try checkFamily(.mips32);
}

test "every core of the loongarch64 family resolves to sane parameters" {
    try checkFamily(.loongarch64);
}

test "every core of the sparc family resolves to sane parameters" {
    try checkFamily(.sparc);
}

test "every core of the wasm family resolves to sane parameters" {
    try checkFamily(.wasm);
}

test "every core of the generic family resolves to sane parameters" {
    try checkFamily(.generic);
}

test "generic follows the GENERIC rule for a 512 KB L2 on x86" {
    var vt: Vtable = .{};
    installGeneric(&vt, .x86, native_buffer_size, .{});
    // min(512 >> 6, 16) = 8
    try std.testing.expectEqual(@as(i32, 448), vt.sgemm_p);
    try std.testing.expectEqual(@as(i32, 224), vt.dgemm_p);
    try std.testing.expectEqual(@as(i32, 112), vt.zgemm_p);
    try std.testing.expectEqual(@as(i32, 128), vt.dgemm_q);
    const packed_a: i64 = (448 * 128 * 4 + 0xffff) & ~@as(i64, 0xffff);
    const want: i64 = (@divTrunc(native_buffer_size - packed_a, 128 * 4) - 15) & ~@as(i64, 15);
    try std.testing.expectEqual(@as(i32, @intCast(want)), vt.sgemm_r);
}

test "the block factor scales a computed p in 8s, and the L2 size moves it" {
    var vt: Vtable = .{};
    installGeneric(&vt, .x86, native_buffer_size, .{ .block_factor = 50 });
    try std.testing.expectEqual(@as(i32, 224), vt.sgemm_p);
    installGeneric(&vt, .x86, native_buffer_size, .{ .block_factor = 1 });
    // clamped to 10 percent: 44 in 8s is 40
    try std.testing.expectEqual(@as(i32, 40), vt.sgemm_p);
    installGeneric(&vt, .x86, native_buffer_size, .{ .l2_kb = 128 });
    try std.testing.expectEqual(@as(i32, 112), vt.sgemm_p);
}

test "r derives from the core's own BUFFER_SIZE, a fallback from its family's" {
    const mib: i64 = 1 << 20;
    inline for (.{
        .{ .x86_64, .haswell, 128 },    .{ .x86, .katmai, 16 },      .{ .arm64, .neoversen1, 128 },
        .{ .arm64, .cortexa53, 32 },    .{ .power64, .power9, 256 }, .{ .power64, .ppc440, 256 },
        .{ .mips64, .loongson3r3, 64 }, .{ .s390x, .z14, 128 },      .{ .wasm, .wasm128_generic, 16 },
    }) |c| try std.testing.expectEqual(c[2] * mib, tableFor(c[0], c[1]).?.buffer_size);
    var vt: Vtable = .{};
    installCore(&vt, .x86_64, .haswell, .{});
    const t = tableFor(.x86_64, .haswell).?;
    try std.testing.expectEqual(derivedR(vt.sgemm_p, vt.sgemm_q, 4, t), vt.sgemm_r);
    const packed_a = (@as(i64, vt.sgemm_p) * vt.sgemm_q * 4 + t.offset_a + t.align_) & ~@as(i64, t.align_);
    try std.testing.expectEqual((@divTrunc(128 * mib - packed_a, @as(i64, vt.sgemm_q) * 4) - 15) & ~@as(i64, 15), vt.sgemm_r);
    // `p6` has no table: the generic set on x86's 16 MiB
    installCore(&vt, .x86, .p6, .{});
    var g = data.generic(.generic).?;
    g.buffer_size = 16 * mib;
    try std.testing.expectEqual(derivedR(vt.sgemm_p, vt.sgemm_q, 4, g), vt.sgemm_r);
}

test "the generic table takes the target's BUFFER_SIZE" {
    const mib: i64 = 1 << 20;
    inline for (.{
        .{ .x86_64, 64, 128 },  .{ .x86, 32, 16 },     .{ .aarch64, 64, 32 },  .{ .arm, 32, 32 },
        .{ .powerpc, 32, 256 }, .{ .sparc64, 64, 32 }, .{ .mips64el, 64, 64 }, .{ .mips64el, 32, 32 },
        .{ .wasm32, 32, 16 },   .{ .alpha, 64, 32 },   .{ .csky, 32, 32 },     .{ .m68k, 32, 128 },
    }) |c| try std.testing.expectEqual(c[2] * mib, genericBufferSize(c[0], c[1]));
}

test "the generic table follows parameter.c's L2 rule on x86 only" {
    inline for (.{
        .{ .x86_64, .x86 },    .{ .x86, .x86 },      .{ .aarch64, .other },           .{ .arm, .other },
        .{ .powerpc, .other }, .{ .mips, .mips },    .{ .mipsel, .mips },             .{ .mips64el, .other },
        .{ .riscv64, .other }, .{ .wasm32, .other }, .{ .loongarch64, .loongarch64 }, .{ .loongarch32, .other },
        .{ .sparc64, .other }, .{ .s390x, .other },  .{ .m68k, .other },
    }) |c| try std.testing.expectEqual(@as(GenericBranch, c[1]), genericBranch(c[0]));
    var vt: Vtable = .{};
    const env: Env = .{ .l2_kb = 4096, .block_factor = 50 };
    // x86: 56 * min(4096 >> 6, 16) at half
    installGeneric(&vt, .x86, 16 << 20, env);
    try std.testing.expectEqual(@as(i32, 448), vt.sgemm_p);
    // elsewhere `parameter.c:46-48`'s 128s, whatever the L2, factor and buffer
    installGeneric(&vt, .other, 32 << 20, env);
    inline for ([_]Gemm{ .s, .d, .q, .c, .z, .x }) |g| {
        inline for ([_]Param{ .p, .q, .r }) |f| try std.testing.expectEqual(@as(i32, 128), slots.paramPtr(&vt, g, f).*);
    }
    // MIPS and LoongArch: `param.h:4456-4483`'s sets, which have no `q` or `x`
    const want = [4][3]i32{ .{ 128, 240, 12288 }, .{ 128, 120, 8192 }, .{ 96, 120, 4096 }, .{ 64, 120, 4096 } };
    inline for (.{ .mips, .loongarch64 }) |br| {
        installGeneric(&vt, br, 32 << 20, env);
        inline for ([_]Gemm{ .s, .d, .c, .z }, want) |g, w| {
            try std.testing.expectEqual(w, [3]i32{ slots.paramPtr(&vt, g, .p).*, slots.paramPtr(&vt, g, .q).*, slots.paramPtr(&vt, g, .r).* });
        }
        try std.testing.expectEqual(@as(i32, 0), vt.qgemm_q);
    }
    // a core with no table takes its architecture's branch
    installCore(&vt, .power64, .unknown, env);
    try std.testing.expectEqual(@as(i32, 128), vt.dgemm_r);
    installCore(&vt, .mips32, .unknown, env);
    try std.testing.expectEqual(@as(i32, 8192), vt.dgemm_r);
}

test "fixed sets are taken as they are; x86 derives r, others keep the header's" {
    var vt: Vtable = .{};
    installCore(&vt, .x86_64, .haswell, .{});
    try std.testing.expectEqual(@as(i32, 320), vt.sgemm_p);
    try std.testing.expectEqual(@as(i32, 8), vt.sgemm_unroll_m);
    try std.testing.expect(vt.dgemm_r != 13824);
    installCore(&vt, .mips64, .mips64_generic, .{});
    try std.testing.expectEqual(@as(i32, 128), vt.sgemm_p);
    try std.testing.expectEqual(@as(i32, 12288), vt.sgemm_r);
}

test "mips64 takes setparam-ref.c's dgemm r of 640 on every core; mips32 keeps its own" {
    var vt: Vtable = .{};
    inline for (.{ .mips64_generic, .loongson3r3, .loongson3r4, .i6400, .p6600 }) |c| {
        installCore(&vt, .mips64, c, .{});
        try std.testing.expectEqual(@as(i32, 640), vt.dgemm_r);
    }
    // the header's 8192 and 12288 before; only `d` moves
    installCore(&vt, .mips64, .mips64_generic, .{});
    try std.testing.expectEqual(@as(i32, 12288), vt.sgemm_r);
    try std.testing.expectEqual(@as(i32, 4096), vt.zgemm_r);
    try std.testing.expectEqual(@as(i32, 128), vt.dgemm_p);
    try std.testing.expectEqual(@as(i32, 120), vt.dgemm_q);
    // `ARCH_MIPS` has no such block: `r` from the buffer
    installCore(&vt, .mips32, .loongson3r3, .{});
    try std.testing.expectEqual(derivedR(44, 92, 8, tableFor(.mips32, .loongson3r3).?), vt.dgemm_r);
}

test "x86 3M takes the real precision's p and q where the header has none; xgemm3m always q's" {
    var vt: Vtable = .{};
    // no 3M macros: `setparam-ref.c:1594, 1602, 2109, 2117`
    inline for (.{ .core2, .nehalem }) |c| {
        installCore(&vt, .x86_64, c, .{ .l2_kb = 4096 });
        const t = tableFor(.x86_64, c).?;
        try std.testing.expect(vt.cgemm3m_p > 0);
        try std.testing.expectEqual(vt.sgemm_p, vt.cgemm3m_p);
        try std.testing.expectEqual(vt.sgemm_q, vt.cgemm3m_q);
        try std.testing.expectEqual(derivedR(vt.cgemm3m_p, vt.cgemm3m_q, 8, t), vt.cgemm3m_r);
        try std.testing.expectEqual(vt.dgemm_p, vt.zgemm3m_p);
        try std.testing.expectEqual(vt.dgemm_q, vt.zgemm3m_q);
        try std.testing.expectEqual(derivedR(vt.zgemm3m_p, vt.zgemm3m_q, 16, t), vt.zgemm3m_r);
    }
    // the header's 320 and 320 stand; `x3m` takes `q`'s 504 and 128 over the
    // header's 112 and 224 (`:1609, 2122`)
    installCore(&vt, .x86_64, .haswell, .{});
    try std.testing.expectEqual(@as(i32, 320), vt.cgemm3m_p);
    try std.testing.expectEqual(@as(i32, 320), vt.cgemm3m_q);
    try std.testing.expectEqual(@as(i32, 504), vt.xgemm3m_p);
    try std.testing.expectEqual(@as(i32, 128), vt.xgemm3m_q);
    // the AMD override moves `s`, not the header's 3M
    installCore(&vt, .x86_64, .skylakex, .{ .amd_avx512 = true, .l2_kb = 1024, .l3_kb = 32768 });
    try std.testing.expectEqual(@as(i32, 384), vt.sgemm_p);
    try std.testing.expectEqual(@as(i32, 320), vt.cgemm3m_p);
}

test "3M stays 0 outside x86, as no DYNAMIC_ARCH build there defines USE_GEMM3M" {
    var vt: Vtable = .{};
    inline for (.{
        .{ .arm64, .neoversen1 }, .{ .mips64, .loongson3r4 }, .{ .riscv64, .riscv64_generic },
        .{ .power64, .power9 },   .{ .s390x, .z14 },          .{ .loongarch64, .la464 },
    }) |c| {
        installCore(&vt, c[0], c[1], .{});
        inline for (.{ "c", "z", "x" }) |pre| {
            inline for (.{ "_p", "_q", "_r" }) |f| try std.testing.expectEqual(@as(i32, 0), @field(vt, pre ++ "gemm3m" ++ f));
        }
    }
}

test "only x86 rounds p up to the unroll; Blocking takes whole panels of an unrounded one" {
    var vt: Vtable = .{};
    // the header's 44 on an unroll of 8, which OpenBLAS keeps (`setparam-ref.c:1101`;
    // static mips32 builds take the header's)
    inline for (.{ .mips64, .mips32 }) |f| {
        installCore(&vt, f, .loongson3r4, .{});
        try std.testing.expectEqual(@as(i32, 8), vt.dgemm_unroll_m);
        try std.testing.expectEqual(@as(i32, 44), vt.dgemm_p);
        try std.testing.expectEqual(@as(usize, 40), driver.Blocking.of(f64, &vt).p);
    }
    // core2's rule: 92 * (512 >> 9) + 8 = 100, rounded to 104; `q`'s
    // 23 + 4 = 27 stays (`:2191-2195`)
    installCore(&vt, .x86_64, .core2, .{ .l2_kb = 512 });
    try std.testing.expectEqual(@as(i32, 104), vt.sgemm_p);
    try std.testing.expectEqual(@as(i32, 104), vt.cgemm3m_p);
    try std.testing.expectEqual(@as(i32, 27), vt.qgemm_p);
}

test "half and bfloat16 take the header's r outside x86, and only what each family's block sets" {
    var vt: Vtable = .{};
    // `b` and `sb` only: `setparam-ref.c:1027, 1198, 1337, 1369`
    inline for (.{ .{ .arm64, .neoversen1 }, .{ .loongarch64, .la464 }, .{ .power64, .power10 }, .{ .s390x, .z14 } }) |c| {
        installCore(&vt, c[0], c[1], .{});
        const t = tableFor(c[0], c[1]).?;
        try std.testing.expectEqual(t.sb.r.fixed, vt.sbgemm_r);
        try std.testing.expectEqual(t.b.r.fixed, vt.bgemm_r);
        try std.testing.expectEqual(t.sb.p.fixed, vt.sbgemm_p);
        inline for ([_]Param{ .p, .q, .r }) |f| try std.testing.expectEqual(@as(i32, 0), slots.paramPtr(&vt, .sh, f).*);
    }
    installCore(&vt, .arm64, .neoversen1, .{});
    try std.testing.expectEqual(@as(i32, 256), vt.sbgemm_r);
    // riscv64 also `sh`, whose header `r` is 240 (`:1409`)
    installCore(&vt, .riscv64, .riscv64_generic, .{});
    try std.testing.expectEqual(@as(i32, 240), vt.shgemm_r);
    try std.testing.expectEqual(@as(i32, 256), vt.sbgemm_r);
    try std.testing.expectEqual(@as(i32, 256), vt.bgemm_r);
    // mips64 none
    installCore(&vt, .mips64, .loongson3r4, .{});
    inline for ([_]Gemm{ .sh, .b, .sb }) |g| {
        inline for ([_]Param{ .p, .q, .r }) |f| try std.testing.expectEqual(@as(i32, 0), slots.paramPtr(&vt, g, f).*);
    }
    // x86 derives it (`:2202-2217`)
    installCore(&vt, .x86_64, .haswell, .{});
    try std.testing.expectEqual(derivedR(256, 256, 4, tableFor(.x86_64, .haswell).?), vt.sbgemm_r);
}

test "old x86 cores follow the L2 size" {
    var vt: Vtable = .{};
    installCore(&vt, .x86_64, .core2, .{ .l2_kb = 4096 });
    // 92 * (4096 >> 9) + 8 = 744
    try std.testing.expectEqual(@as(i32, 744), vt.sgemm_p);
    installCore(&vt, .x86_64, .atom, .{ .l2_kb = 4096 });
    try std.testing.expectEqual(@as(i32, 256), vt.sgemm_p);
}

test "old x86 cores take their getarch entry's subarch blocks" {
    var vt: Vtable = .{};
    // `PENTIUM3` with `HAVE_SSE`; 64 * (512 >> 7)
    installCore(&vt, .x86, .katmai, .{ .l2_kb = 512 });
    try std.testing.expectEqual(@as(i32, 8), vt.sgemm_unroll_m);
    try std.testing.expectEqual(@as(i32, 4), vt.cgemm_unroll_m);
    try std.testing.expectEqual(@as(i32, 256), vt.sgemm_p);
    try std.testing.expectEqual(@as(i32, 256), vt.dgemm_q);
    // `PENTIUM2` (`FORCE_P2`), no SSE; static `p` stays 128
    installCore(&vt, .x86, .p5, .{});
    try std.testing.expectEqual(@as(i32, 4), vt.sgemm_unroll_m);
    try std.testing.expectEqual(@as(i32, 2), vt.cgemm_unroll_m);
    try std.testing.expectEqual(@as(i32, 128), vt.sgemm_p);
    // `PENTIUMM` with `CORE_YONAH`, on Katmai's L2 rule; 16 * (512 >> 7)
    installCore(&vt, .x86, .yonah, .{});
    try std.testing.expectEqual(@as(i32, 4), vt.dgemm_unroll_n);
    try std.testing.expectEqual(@as(i32, 2), vt.zgemm_unroll_n);
    try std.testing.expectEqual(@as(i32, 64), vt.zgemm_p);
    // Athlon's rule goes over its header's fixed 208; 64 * (1024 >> 7)
    installCore(&vt, .x86, .athlon, .{ .l2_kb = 1024 });
    try std.testing.expectEqual(@as(i32, 512), vt.sgemm_p);
    try std.testing.expectEqual(@as(i32, 256), vt.dgemm_p);
    // `CORE_PRESCOTT`'s offsets follow `__64BIT__`
    installCore(&vt, .x86, .prescott, .{});
    try std.testing.expectEqual(@as(i32, 128), vt.offsetA);
    installCore(&vt, .x86_64, .prescott, .{});
    try std.testing.expectEqual(@as(i32, 256), vt.offsetB);
}

test "32-bit builds take Makefile.system's fallback entries" {
    const pairs = .{
        .{ .x86, .haswell, .nehalem },     .{ .x86, .skylakex, .nehalem },  .{ .x86, .sapphirerapids, .nehalem },
        .{ .x86, .sandybridge, .nehalem }, .{ .x86, .zen, .barcelona },     .{ .x86, .bulldozer, .barcelona },
        .{ .x86, .excavator, .barcelona }, .{ .power32, .power8, .power6 }, .{ .power32, .power10, .power6 },
        .{ .power32le, .power9, .power6 },
    };
    inline for (pairs) |p| {
        try std.testing.expectEqual(tableFor(p[0], p[2]), tableFor(p[0], p[1]));
        try std.testing.expectEqual(dtbFor(p[0], p[2]), dtbFor(p[0], p[1]));
    }
    var vt: Vtable = .{};
    installCore(&vt, .x86, .haswell, .{});
    try std.testing.expectEqual(@as(i32, 504), vt.sgemm_p);
    try std.testing.expectEqual(@as(i32, 512), vt.sgemm_q);
    try std.testing.expectEqual(@as(i32, 256), vt.dgemm_q);
    installCore(&vt, .x86, .zen, .{});
    try std.testing.expectEqual(@as(i32, 448), vt.sgemm_p);
    try std.testing.expectEqual(@as(i32, 224), vt.sgemm_q);
    // 64-bit keeps the cores' own
    installCore(&vt, .x86_64, .haswell, .{});
    try std.testing.expectEqual(@as(i32, 320), vt.sgemm_p);
    try std.testing.expect(!std.meta.eql(tableFor(.power64, .power8), tableFor(.power64, .power6)));
}

/// The 32-bit x86 objects, as a `select.served` context.
const X86Objects = struct {
    pub const arch: std.Target.Cpu.Arch = .x86;
    pub const fallback: ?coretypes.X86CoreType = null;
    pub fn ok(_: X86Objects, c: coretypes.X86CoreType) bool {
        return std.mem.findScalar(coretypes.X86CoreType, corespec.kernelCores(.x86), c) != null;
    }
};

test "a 32-bit x86 core's serving object has the core's blocking" {
    for (std.enums.values(coretypes.X86CoreType)) |c| {
        if (!select.selectable(c)) continue;
        const k = select.served(coretypes.X86CoreType, c, X86Objects{}) orelse continue;
        try std.testing.expectEqual(tableFor(.x86, c), tableFor(.x86, k));
        try std.testing.expectEqual(dtbFor(.x86, c), dtbFor(.x86, k));
    }
}

test "32-bit x86 serves the AVX cores and Zen from their fallbacks' objects" {
    const X = coretypes.X86CoreType;
    const pairs = [_][2]X{
        .{ .sandybridge, .nehalem },  .{ .haswell, .nehalem },        .{ .skylakex, .nehalem },
        .{ .cooperlake, .nehalem },   .{ .sapphirerapids, .nehalem }, .{ .bulldozer, .barcelona },
        .{ .piledriver, .barcelona }, .{ .steamroller, .barcelona },  .{ .excavator, .barcelona },
        .{ .zen, .barcelona },
    };
    for (pairs) |p| {
        try std.testing.expect(std.mem.findScalar(X, corespec.kernelCores(.x86), p[0]) == null);
        try std.testing.expectEqual(@as(?X, p[1]), select.served(X, p[0], X86Objects{}));
        try std.testing.expectEqual(tableFor(.x86, p[0]), tableFor(.x86, p[1]));
        try std.testing.expectEqual(dtbFor(.x86, p[0]), dtbFor(.x86, p[1]));
        // 64-bit builds keep them
        try std.testing.expect(std.mem.findScalar(X, corespec.kernelCores(.x86_64), p[0]) != null);
    }
}

test "PPC970's unrolls follow the byte order" {
    var vt: Vtable = .{};
    inline for (.{ .{ .power64, 16, 8 }, .{ .power64be, 4, 2 }, .{ .power32, 4, 2 }, .{ .power32le, 16, 8 } }) |f| {
        installCore(&vt, f[0], .ppc970, .{});
        try std.testing.expectEqual(@as(i32, f[1]), vt.sgemm_unroll_m);
        try std.testing.expectEqual(@as(i32, f[2]), vt.cgemm_unroll_m);
        try std.testing.expectEqual(@as(i32, 4), vt.dgemm_unroll_m);
    }
}

test "LoongArch's old names are their getarch aliases" {
    var vt: Vtable = .{};
    installCore(&vt, .loongarch64, .loongson3r5, .{ .threads = 1, .l3_kb = 32 * 1024 });
    try std.testing.expectEqual(@as(i32, 8192), vt.sgemm_r);
    try std.testing.expectEqual(@as(i32, 16), vt.sgemm_unroll_m);
}

test "LA464 picks its set by threads and L3 size" {
    var vt: Vtable = .{};
    installCore(&vt, .loongarch64, .la464, .{ .threads = 1, .l3_kb = 32 * 1024 });
    try std.testing.expectEqual(@as(i32, 8192), vt.sgemm_r);
    try std.testing.expectEqual(@as(i32, 289), vt.dgemm_q);
    installCore(&vt, .loongarch64, .la464, .{ .threads = 4, .l3_kb = 16 * 1024 });
    try std.testing.expectEqual(@as(i32, 2048), vt.sgemm_r);
    try std.testing.expectEqual(@as(i32, 738), vt.dgemm_r);
}

test "detection gives usable cache sizes on this machine" {
    const env = detect();
    try std.testing.expect(env.l2_kb > 0);
    try std.testing.expect(env.l3_kb >= 0);
}
