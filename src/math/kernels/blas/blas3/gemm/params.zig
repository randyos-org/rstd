//! The gemm blocking parameters, put in the vtable as a core's kernels are
//! installed, as `kernel/setparam-ref.c` fills OpenBLAS's table: `param.h`'s
//! constants (`params_data/`) and the rules it computes (`notes.md`).

const std = @import("std");
const Vtable = @import("../../Vtable.zig");
const Config = @import("../../Config.zig");
const cache = @import("../../cpu/detect/cache.zig");
const slots = @import("../../slots.zig");
const plook = @import("params/lookup.zig");

const Gemm = slots.Gemm;
pub const CoreOf = plook.CoreOf;
pub const Family = plook.Family;
pub const GenericBranch = plook.GenericBranch;
pub const dtbFor = plook.dtbFor;
pub const genericBranch = plook.genericBranch;
pub const genericBufferSize = plook.genericBufferSize;
pub const installCore = plook.installCore;
pub const installGeneric = plook.installGeneric;
pub const native_branch = plook.native_branch;
pub const native_buffer_size = plook.native_buffer_size;
pub const tableFor = plook.tableFor;
pub const tables = plook.tables;

/// How `p`, `q` or `r` comes about.
pub const Size = union(enum) {
    /// A constant of `param.h`.
    fixed: i32,
    /// The header names a run-time variable, or has no macro; see `Rule`.
    dynamic,
    /// `r` derived from the packing buffer once `p` and `q` are known; also
    /// an `r` the header has no macro for.
    from_buffer,
};

/// One precision's parameters. A shape the generator leaves out (no `P` or
/// unroll macro) defaults to what it writes for a missing macro.
pub const Shape = struct {
    unroll_m: i32 = 0,
    unroll_n: i32 = 0,
    /// `GEMM_UNROLL_MN`.
    unroll_mn: i32 = 0,
    /// `BGEMM_ALIGN_K` and `SBGEMM_ALIGN_K`.
    align_k: i32 = 0,
    /// `<P>GEMM_DEFAULT_P_BASE`, where the header has one.
    p_base: i32 = 0,
    p: Size = .dynamic,
    q: Size = .dynamic,
    r: Size = .from_buffer,
};

/// A core's whole table, one shape per precision: the four the BLAS names, the
/// extended ones, half and bfloat16, and the 3M complex variants.
pub const Table = struct {
    s: Shape = .{},
    d: Shape = .{},
    q: Shape = .{},
    c: Shape = .{},
    z: Shape = .{},
    x: Shape = .{},
    sh: Shape = .{},
    b: Shape = .{},
    sb: Shape = .{},
    c3m: Shape = .{},
    z3m: Shape = .{},
    x3m: Shape = .{},
    offset_a: i32 = 0,
    offset_b: i32 = 0,
    /// `GEMM_DEFAULT_ALIGN`.
    align_: i32 = 0xffff,
    switch_ratio: i32 = 2,
    divide_rate: i32 = 2,
    divide_limit: i32 = 0,
    preferred_size: i32 = 1,
    /// `SYMV_P`: the block of the generic symv and hemv.
    symv_p: i32 = 0,
    /// `BUFFER_SIZE`, the packing buffer `r` derives from, as the core's
    /// `DYNAMIC_ARCH` build sees `common_<arch>.h`: 128 MiB on x86_64
    /// (`common_x86_64.h:260`), s390x (`common_zarch.h:129`) and arm64's
    /// `NEOVERSEN1`, `N2`, `V1`, `A64FX`, `ARMV8SVE`, `ARMV9SME`
    /// (`common_arm64.h:180-181`); 32 MiB on the other arm64 cores (`:183`), arm
    /// (`common_arm.h:141`), riscv64 (`common_riscv64.h:89`), mips32
    /// (`common_mips.h:105`), LoongArch (`common_loongarch64.h:315`) and SPARC
    /// (`common_sparc.h:240`); 64 MiB on mips64 (`common_mips64.h:254`); 256 MiB
    /// on every POWER core (`common_power.h:923-926`, over the static 2 to 256
    /// MiB of `:914-922`); 16 MiB on x86 (`common_x86.h:207`) and wasm
    /// (`common_wasm.h:87`). The generic table takes the target's
    /// (`native_buffer_size`): 32 MiB on Alpha (`common_alpha.h:128`) and C-SKY
    /// (`common_csky.h:53`).
    buffer_size: i64,
};

/// `p = min(cap, l2_kb >> shift) * mul + add` for `s, d, c, z, q, x` (`cap` 0 is
/// no cap; a zero `mul` with an `add` is a constant).
pub const Affine = struct {
    mul: [6]i32,
    shift: u5 = 7,
    add: [6]i32 = .{ 0, 0, 0, 0, 0, 0 },
    cap: i32 = 0,
};

/// `setparam-ref.c`'s rule for a core whose `p` the header leaves open.
pub const Rule = union(enum) {
    none,
    affine: Affine,
    la464,
    /// `LA264`: `z` by processor id.
    la264,
};

/// The machine the rules read: cache sizes and thread count.
pub const Env = struct {
    /// KB.
    l2_kb: i32 = 512,
    l3_kb: i32 = 0,
    /// `OPENBLAS_NUM_THREADS`; 1 is single-threaded.
    threads: i32 = 0,
    /// x86 only: an AMD part with AVX-512 (Zen 4 and later).
    amd_avx512: bool = false,
    /// LoongArch processor id (`cpucfg` 0).
    prid: u32 = 0,
    /// `OPENBLAS_BLOCK_FACTOR`, percent; 0 leaves `p` alone.
    block_factor: i32 = 0,

    pub fn fromConfig(config: Config) Env {
        var env = detect();
        if (config.l2_size_override > 0) env.l2_kb = config.l2_size_override;
        env.threads = config.num_threads;
        env.block_factor = config.block_factor;
        return env;
    }
};

/// This machine's cache sizes (`cpu/detect/cache.zig`), the rest at their
/// defaults.
pub fn detect() Env {
    const c = cache.detect();
    return .{ .l2_kb = c.l2_kb, .l3_kb = c.l3_kb, .amd_avx512 = c.amd_avx512, .prid = c.prid };
}

// ---- the rules -------------------------------------------------------------

/// `setparam-ref.c`'s L2 rules of the x86 cores whose `p` is open, the LA464
/// override, and `parameter.c`'s `GENERIC` rule (`:451`), which is x86's
/// (`installGeneric` takes `.none` off x86).
pub fn ruleFor(comptime family: Family, comptime core: CoreOf(family)) Rule {
    return switch (family) {
        .generic => .{ .affine = .{ .mul = .{ 56, 28, 28, 14, 14, 7 }, .shift = 6, .cap = 16 } },
        // the old names are `getarch.c` aliases of these
        .loongarch64 => switch (core) {
            .la464, .loongson3r5 => .la464,
            .la264, .loongson2k1000 => .la264,
            else => .none,
        },
        .x86, .x86_64 => switch (core) {
            // `setparam-ref.c:1612`, over Athlon's fixed header `p`
            .katmai, .coppermine, .banias, .yonah, .athlon => .{ .affine = .{ .mul = .{ 64, 32, 32, 16, 16, 8 } } },
            // `setparam-ref.c` has no block for it: `parameter.c:46`'s `DEFAULT_GEMM_P`
            .p5 => .{ .affine = .{ .mul = .{ 0, 0, 0, 0, 0, 0 }, .add = .{ 128, 128, 128, 128, 128, 128 } } },
            .northwood => .{ .affine = .{ .mul = .{ 96, 48, 48, 24, 24, 12 } } },
            .atom => .{ .affine = .{ .mul = .{ 0, 0, 0, 0, 0, 0 }, .add = .{ 256, 128, 128, 64, 64, 32 } } },
            .prescott => .{ .affine = .{ .mul = .{ 56, 28, 28, 14, 14, 7 } } },
            .core2 => .{ .affine = .{ .mul = .{ 92, 46, 46, 23, 23, 11 }, .shift = 9, .add = .{ 8, 8, 4, 4, 4, 2 } } },
            .penryn, .dunnington => .{ .affine = .{ .mul = .{ 42, 42, 21, 21, 42, 21 }, .shift = 9, .add = .{ 8, 8, 4, 4, 8, 4 } } },
            .opteron, .opteron_sse3 => .{ .affine = .{ .mul = .{ 56, 28, 28, 14, 14, 7 }, .add = .{ 224, 112, 112, 56, 56, 28 } } },
            else => .none,
        },
        else => .none,
    };
}

/// One precision's set of `p`, `q`, `r`.
const Resolved = struct { p: i32, q: i32, r: i32 };

/// `LA464`'s sets, `[single thread][L3 is 32 MB][s, d, c, z]` (3C5000/3D5000 have
/// 32 MB).
const la464_sets = [2][2][4]Resolved{
    // multi thread
    .{
        .{ .{ .p = 256, .q = 384, .r = 2048 }, .{ .p = 112, .q = 300, .r = 738 }, .{ .p = 128, .q = 256, .r = 1024 }, .{ .p = 128, .q = 128, .r = 1024 } },
        .{ .{ .p = 256, .q = 384, .r = 1024 }, .{ .p = 112, .q = 289, .r = 353 }, .{ .p = 128, .q = 256, .r = 512 }, .{ .p = 128, .q = 128, .r = 512 } },
    },
    // single thread
    .{
        .{ .{ .p = 256, .q = 384, .r = 4096 }, .{ .p = 112, .q = 300, .r = 3024 }, .{ .p = 128, .q = 256, .r = 2048 }, .{ .p = 128, .q = 128, .r = 1024 } },
        .{ .{ .p = 256, .q = 384, .r = 8192 }, .{ .p = 112, .q = 289, .r = 4096 }, .{ .p = 128, .q = 256, .r = 4096 }, .{ .p = 128, .q = 128, .r = 2048 } },
    },
};

fn roundUp(v: i32, to: i32) i32 {
    if (to <= 0) return v;
    return @divTrunc(v + to - 1, to) * to;
}

/// Applies `OPENBLAS_BLOCK_FACTOR` the way `parameter.c` does.
fn scaled(p: i32, env: Env) i32 {
    if (env.block_factor <= 0) return p;
    const factor: f64 = @floatFromInt(std.math.clamp(env.block_factor, 10, 200));
    const v: i64 = @trunc(@as(f64, @floatFromInt(p)) * factor * 1.0e-2);
    return @intCast(v & ~@as(i64, 7));
}

fn fixedOf(s: Size) ?i32 {
    return switch (s) {
        .fixed => |v| v,
        else => null,
    };
}

/// `p` of one of the primary precisions (`idx`: s, d, c, z, q, x), not yet
/// rounded; a rule goes over the header, as `setparam-ref.c`'s core blocks do.
fn primaryP(shape: Shape, rule: Rule, idx: usize, env: Env) i32 {
    return switch (rule) {
        .affine => |a| blk: {
            const size = if (a.cap > 0) @min(env.l2_kb >> a.shift, a.cap) else env.l2_kb >> a.shift;
            break :blk scaled(a.mul[idx] * size + a.add[idx], env);
        },
        else => fixedOf(shape.p) orelse 0,
    };
}

/// `r` from the buffer: `setparam-ref.c`'s formula.
pub fn derivedR(p: i32, q: i32, elem: i64, table: Table) i32 {
    if (q <= 0) return 0;
    const align_: i64 = table.align_;
    const packed_a = (@as(i64, p) * q * elem + table.offset_a + align_) & ~align_;
    const r = @divTrunc(table.buffer_size - packed_a, @as(i64, q) * elem) - 15;
    return @intCast(r & ~@as(i64, 15));
}

/// Stores the values of `g`, a `Table` field's gemm, in its vtable fields.
fn setPrec(vt: *Vtable, comptime g: Gemm, sh: Shape, r: Resolved) void {
    slots.paramPtr(vt, g, .p).* = r.p;
    slots.paramPtr(vt, g, .q).* = r.q;
    slots.paramPtr(vt, g, .r).* = r.r;
    slots.paramPtr(vt, g, .unroll_m).* = sh.unroll_m;
    slots.paramPtr(vt, g, .unroll_n).* = sh.unroll_n;
    slots.paramPtr(vt, g, .unroll_mn).* = if (sh.unroll_mn > 0) sh.unroll_mn else @max(sh.unroll_m, sh.unroll_n);
    if (comptime g == .b or g == .sb) slots.paramPtr(vt, g, .align_k).* = sh.align_k;
}

/// Fills the vtable's blocking parameters for a core of `family`.
pub fn install(vt: *Vtable, comptime family: Family, comptime core: CoreOf(family), table: Table, env: Env) void {
    installRule(vt, family, core, table, env, comptime ruleFor(family, core));
}

/// `install` under `rule`.
pub fn installRule(vt: *Vtable, comptime family: Family, comptime core: CoreOf(family), table: Table, env: Env, comptime rule: Rule) void {
    const x86 = comptime (family == .x86 or family == .x86_64);

    vt.offsetA = table.offset_a;
    vt.offsetB = table.offset_b;
    vt.@"align" = table.align_;
    vt.switch_ratio = table.switch_ratio;
    vt.divide_rate = table.divide_rate;
    vt.divide_limit = table.divide_limit;
    vt.preferred_size = table.preferred_size;
    vt.symv_p = table.symv_p;

    // the primaries, whose `p` the rules and overrides decide
    const primaries = [_]Gemm{ .s, .d, .c, .z, .q, .x };
    const elems = [6]i64{ 4, 8, 8, 16, 16, 32 };
    var res: [6]Resolved = undefined;
    inline for (primaries, 0..) |g, i| {
        const sh = @field(table, @tagName(g));
        res[i] = .{ .p = primaryP(sh, rule, i, env), .q = fixedOf(sh.q) orelse 0, .r = 0 };
    }

    // 3M, which only x86 sets (`setparam-ref.c:1590-1609, 2105-2122`): the
    // header's `p` and `q`, else the real precision's `p` before the AMD
    // override and rounding, and its header `q`; `x3m` always `q`'s. The
    // arm64 and mips64 blocks set it under `USE_GEMM3M`, which no
    // `DYNAMIC_ARCH` build defines there (`kernel/Makefile.L3:28-45`;
    // `cmake/arch.cmake:132-135` drops mips64's), so 3M stays 0 elsewhere.
    const three_m = [_]Gemm{ .c3m, .z3m, .x3m };
    var res3m: [3]Resolved = @splat(.{ .p = 0, .q = 0, .r = 0 });
    if (x86) {
        inline for (.{ Gemm.c3m, Gemm.z3m }, 0..) |g, i| {
            const sh = @field(table, @tagName(g));
            res3m[i] = .{ .p = fixedOf(sh.p) orelse res[i].p, .q = fixedOf(sh.q) orelse res[i].q, .r = 0 };
        }
        res3m[2] = .{ .p = res[4].p, .q = res[4].q, .r = 0 };
    }

    if (comptime family == .x86_64) {
        // AMD parts with AVX-512 run the Skylake-X kernels on other blocking
        const skx = comptime switch (core) {
            .skylakex, .cooperlake, .sapphirerapids => true,
            else => false,
        };
        if (skx and env.amd_avx512 and env.l3_kb > 0 and @rem(env.l3_kb, 32768) == 0 and env.l2_kb == 1024) {
            const pq = [4][2]i32{ .{ 384, 512 }, .{ 512, 512 }, .{ 160, 480 }, .{ 176, 256 } };
            inline for (0..4) |i| res[i] = .{ .p = pq[i][0], .q = pq[i][1], .r = 0 };
        }
    }

    // only x86 rounds `p` up to the unroll (`setparam-ref.c:2162-2187`): `s`,
    // `d`, `c`, `z`, and 3M to its own unroll, else the real one's; `q`, `x`
    // and `x3m` only under `QUAD_PRECISION` (`:2191-2195`), off by default
    // (`Makefile.rule:214`)
    if (x86) {
        inline for (.{ Gemm.s, Gemm.d, Gemm.c, Gemm.z }, 0..) |g, i| res[i].p = roundUp(res[i].p, @field(table, @tagName(g)).unroll_m);
        const unroll_of = [2]i32{ table.s.unroll_m, table.d.unroll_m };
        inline for (.{ Gemm.c3m, Gemm.z3m }, 0..) |g, i| {
            const um = @field(table, @tagName(g)).unroll_m;
            res3m[i].p = roundUp(res3m[i].p, if (um > 0) um else unroll_of[i]);
        }
    }

    // x86 always derives `r`, the other families take the header's constant
    inline for (primaries, 0..) |g, i| {
        const fixed = if (x86) null else fixedOf(@field(table, @tagName(g)).r);
        res[i].r = fixed orelse derivedR(res[i].p, res[i].q, elems[i], table);
    }
    if (x86) {
        const elems3m = [3]i64{ 8, 16, 32 };
        inline for (0..3) |i| res3m[i].r = derivedR(res3m[i].p, res3m[i].q, elems3m[i], table);
    }

    // `setparam-ref.c:1111` (`ARCH_MIPS64`), over the header's 8192 or the
    // Loongsons' `dgemm_r`; `s`, `c` and `z` keep the header's
    if (comptime family == .mips64) res[1].r = 640;

    if (comptime std.meta.activeTag(rule) == .la464) {
        const set = la464_sets[@intFromBool(env.threads == 1)][@intFromBool(@divTrunc(env.l3_kb, 1024) == 32)];
        res[0..4].* = set;
    }
    if (comptime std.meta.activeTag(rule) == .la264) {
        // the 2K3000 (processor id 0x0014b020) has its own `z` set
        const set: Resolved = if (env.prid == 0x0014b020) .{ .p = 128, .q = 176, .r = 360 } else .{ .p = 64, .q = 120, .r = 4096 };
        res[3] = set;
    }
    inline for (primaries, 0..) |g, i| setPrec(vt, g, @field(table, @tagName(g)), res[i]);
    inline for (three_m, 0..) |g, i| setPrec(vt, g, @field(table, @tagName(g)), res3m[i]);

    // half and bfloat16: the headers' `p` and `q`; x86 derives `r`
    // (`setparam-ref.c:1568-1575, 2202-2217`), the others take the header's.
    // arm64, LoongArch, POWER and zarch set only `b` and `sb` (`:993-1027,
    // 1193-1198, 1319, 1328-1347, 1360-1379`), riscv64 also `sh`
    // (`:1393-1422`), mips64 none; the static families' are the header's
    // (`common_param.h:1480-1504`).
    inline for (.{ Gemm.sh, Gemm.b, Gemm.sb }) |g| {
        const set = comptime switch (family) {
            .mips64 => false,
            .arm64, .loongarch64, .power64, .power64be, .power32, .power32le, .s390x => g != .sh,
            .x86_64, .x86, .riscv64, .arm, .mips32, .sparc, .wasm, .generic => true,
        };
        const sh = @field(table, @tagName(g));
        const p = fixedOf(sh.p) orelse 0;
        const q = fixedOf(sh.q) orelse 0;
        const fixed = if (x86) null else fixedOf(sh.r);
        setPrec(vt, g, sh, if (set) .{ .p = p, .q = q, .r = fixed orelse derivedR(p, q, 4, table) } else .{ .p = 0, .q = 0, .r = 0 });
    }
}
