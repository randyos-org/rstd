//! The gemm blocking parameters, put in the vtable when a core's kernels are
//! installed, the way OpenBLAS fills its per-core parameter table
//! (`kernel/setparam-ref.c`, the `DYNAMIC_ARCH` builder, which is the model for
//! choosing the core at run time).
//!
//! Per precision the table (`params_data.zig`, generated from `param.h` by
//! `tools/blas_asm/params.zig`) has the unroll sizes, which the kernels fix, and
//! `p` (rows of A packed at a time), `q` (the depth) and `r` (columns of B). Most
//! are constants of `param.h`. What `setparam-ref.c` computes is here:
//!
//!   * `p` of the older x86 cores follows the L2 size (`Affine`), over the header;
//!   * on x86, `r` is always derived from the packing buffer;
//!   * on mips64, `dgemm`'s `r` is 640 for every core;
//!   * half and bfloat16 are set per family, `r` the header's outside x86;
//!   * on x86, `p` of `s`, `d`, `c`, `z` and 3M is rounded up to the unroll;
//!   * on x86 only, the 3M variants take `p` and `q` of the real ones unless
//!     the header says;
//!   * LA464 and the AMD AVX-512 parts override sets by cache size and threads.
//!
//! `OPENBLAS_BLOCK_FACTOR` (`Config.block_factor`) scales a computed `p` as
//! `driver/others/parameter.c` does there; OpenBLAS ignores it in dynamic builds,
//! here it applies wherever `p` comes from a rule.

const std = @import("std");
const sysinfo = @import("builtin");
const Vtable = @import("../../Vtable.zig");
const Config = @import("../../Config.zig");
const data = @import("params_data.zig");
const coretypes = @import("../../cpu/_coretypes.zig");
const corespec = @import("../../cpu/_corespec.zig");
const select = @import("../../cpu/_select.zig");
const driver = @import("../driver.zig");

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
    /// (`common_wasm.h:87`).
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

/// The cache sizes of this machine where an instruction gives them: cpuid on
/// x86 (`get_l2_size` / `get_l3_size`) and cpucfg on LoongArch. Elsewhere, and
/// where the answer is zero, the defaults stand (OpenBLAS assumes 256 KB on x86
/// and 512 KB on RISC-V, which reads sysfs and so needs an OS).
pub fn detect() Env {
    var env: Env = .{};
    switch (sysinfo.target.cpu.arch) {
        .x86_64, .x86 => {
            const id = cpuid(0x80000006, 0);
            const l2: i32 = @intCast((id.ecx >> 16) & 0xffff);
            env.l2_kb = if (l2 > 0) l2 else 256;
            env.l3_kb = @intCast(((id.edx >> 18) & 0x3fff) * 512);
            const vendor = cpuid(0, 0);
            const leaf7 = cpuid(7, 0);
            // "Auth" of "AuthenticAMD"
            env.amd_avx512 = vendor.ebx == 0x68747541 and (leaf7.ebx & (1 << 16)) != 0;
        },
        .loongarch64 => {
            const cfg = cpucfg(0x14);
            const mb = ((cfg & 0xffff) + 1) * (@as(u64, 1) << @intCast((cfg >> 16) & 0xff)) * (@as(u64, 1) << @intCast((cfg >> 24) & 0x7f)) / 1024 / 1024;
            env.l3_kb = @intCast(mb * 1024);
            env.prid = @truncate(cpucfg(0));
        },
        else => {},
    }
    return env;
}

const CpuidResult = struct { eax: u32, ebx: u32, ecx: u32, edx: u32 };

fn cpuid(leaf: u32, sub: u32) CpuidResult {
    switch (sysinfo.target.cpu.arch) {
        .x86_64, .x86 => {
            var a: u32 = undefined;
            var b: u32 = undefined;
            var c: u32 = undefined;
            var d: u32 = undefined;
            asm volatile ("cpuid"
                : [a] "={eax}" (a),
                  [b] "={ebx}" (b),
                  [c] "={ecx}" (c),
                  [d] "={edx}" (d),
                : [leaf] "{eax}" (leaf),
                  [sub] "{ecx}" (sub),
            );
            return .{ .eax = a, .ebx = b, .ecx = c, .edx = d };
        },
        else => return .{ .eax = 0, .ebx = 0, .ecx = 0, .edx = 0 },
    }
}

fn cpucfg(word: u32) u64 {
    switch (sysinfo.target.cpu.arch) {
        .loongarch64 => return asm volatile ("cpucfg %[ret], %[id]"
            : [ret] "=r" (-> u64),
            : [id] "r" (word),
        ),
        else => return 0,
    }
}

// ---- the rules -------------------------------------------------------------

/// `setparam-ref.c`'s L2 rules of the x86 cores whose `p` is open, the LA464
/// override, and the `parameter.c` `GENERIC` rule.
pub fn ruleFor(comptime family: Family, comptime core_name: []const u8) Rule {
    const eq = std.mem.eql;
    if (comptime family == .generic) return .{ .affine = .{ .mul = .{ 56, 28, 28, 14, 14, 7 }, .shift = 6, .cap = 16 } };
    // the old names are `getarch.c` aliases of these
    if (comptime family == .loongarch64 and (eq(u8, core_name, "la464") or eq(u8, core_name, "loongson3r5"))) return .la464;
    if (comptime family == .loongarch64 and (eq(u8, core_name, "la264") or eq(u8, core_name, "loongson2k1000"))) return .la264;
    if (comptime family != .x86 and family != .x86_64) return .none;
    // `setparam-ref.c:1600`, over Athlon's fixed header `p`
    if (comptime eq(u8, core_name, "katmai") or eq(u8, core_name, "coppermine") or eq(u8, core_name, "banias") or
        eq(u8, core_name, "yonah") or eq(u8, core_name, "athlon"))
        return .{ .affine = .{ .mul = .{ 64, 32, 32, 16, 16, 8 } } };
    // `setparam-ref.c` has no block for it: `parameter.c:46`'s `DEFAULT_GEMM_P`
    if (comptime eq(u8, core_name, "p5"))
        return .{ .affine = .{ .mul = .{ 0, 0, 0, 0, 0, 0 }, .add = .{ 128, 128, 128, 128, 128, 128 } } };
    if (comptime eq(u8, core_name, "northwood")) return .{ .affine = .{ .mul = .{ 96, 48, 48, 24, 24, 12 } } };
    if (comptime eq(u8, core_name, "atom")) return .{ .affine = .{ .mul = .{ 0, 0, 0, 0, 0, 0 }, .add = .{ 256, 128, 128, 64, 64, 32 } } };
    if (comptime eq(u8, core_name, "prescott")) return .{ .affine = .{ .mul = .{ 56, 28, 28, 14, 14, 7 } } };
    if (comptime eq(u8, core_name, "core2"))
        return .{ .affine = .{ .mul = .{ 92, 46, 46, 23, 23, 11 }, .shift = 9, .add = .{ 8, 8, 4, 4, 4, 2 } } };
    if (comptime eq(u8, core_name, "penryn") or eq(u8, core_name, "dunnington"))
        return .{ .affine = .{ .mul = .{ 42, 42, 21, 21, 42, 21 }, .shift = 9, .add = .{ 8, 8, 4, 4, 8, 4 } } };
    if (comptime eq(u8, core_name, "opteron") or eq(u8, core_name, "opteron_sse3"))
        return .{ .affine = .{ .mul = .{ 56, 28, 28, 14, 14, 7 }, .add = .{ 224, 112, 112, 56, 56, 28 } } };
    return .none;
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
fn derivedR(p: i32, q: i32, elem: i64, table: Table) i32 {
    if (q <= 0) return 0;
    const align_: i64 = table.align_;
    const packed_a = (@as(i64, p) * q * elem + table.offset_a + align_) & ~align_;
    const r = @divTrunc(table.buffer_size - packed_a, @as(i64, q) * elem) - 15;
    return @intCast(r & ~@as(i64, 15));
}

/// Stores one precision's values under the vtable's `<name>gemm_*` fields.
fn setPrec(vt: *Vtable, comptime name: []const u8, sh: Shape, r: Resolved) void {
    const prefix = comptime if (std.mem.eql(u8, name, "c3m")) "cgemm3m" else if (std.mem.eql(u8, name, "z3m")) "zgemm3m" else if (std.mem.eql(u8, name, "x3m")) "xgemm3m" else name ++ "gemm";
    @field(vt, prefix ++ "_p") = r.p;
    @field(vt, prefix ++ "_q") = r.q;
    @field(vt, prefix ++ "_r") = r.r;
    @field(vt, prefix ++ "_unroll_m") = sh.unroll_m;
    @field(vt, prefix ++ "_unroll_n") = sh.unroll_n;
    @field(vt, prefix ++ "_unroll_mn") = if (sh.unroll_mn > 0) sh.unroll_mn else @max(sh.unroll_m, sh.unroll_n);
    if (comptime std.mem.eql(u8, name, "b") or std.mem.eql(u8, name, "sb")) @field(vt, prefix ++ "_align_k") = sh.align_k;
}

/// Fills the vtable's blocking parameters for a core of `family`.
pub fn install(vt: *Vtable, comptime family: Family, comptime core_name: []const u8, table: Table, env: Env) void {
    const rule = comptime ruleFor(family, core_name);
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
    const names = .{ "s", "d", "c", "z", "q", "x" };
    const elems = [6]i64{ 4, 8, 8, 16, 16, 32 };
    var res: [6]Resolved = undefined;
    inline for (names, 0..) |name, i| {
        const sh = @field(table, name);
        res[i] = .{ .p = primaryP(sh, rule, i, env), .q = fixedOf(sh.q) orelse 0, .r = 0 };
    }

    // 3M, which only x86 sets (`setparam-ref.c:1578-1597, 2093-2110`): the
    // header's `p` and `q`, else the real precision's `p` before the AMD
    // override and rounding, and its header `q`; `x3m` always `q`'s. The
    // arm64 and mips64 blocks set it under `USE_GEMM3M`, which no
    // `DYNAMIC_ARCH` build defines there (`kernel/Makefile.L3:28-45`;
    // `cmake/arch.cmake:132-135` drops mips64's), so 3M stays 0 elsewhere.
    const names3m = .{ "c3m", "z3m", "x3m" };
    var res3m: [3]Resolved = @splat(.{ .p = 0, .q = 0, .r = 0 });
    if (x86) {
        inline for (.{ "c3m", "z3m" }, 0..) |name, i| {
            const sh = @field(table, name);
            res3m[i] = .{ .p = fixedOf(sh.p) orelse res[i].p, .q = fixedOf(sh.q) orelse res[i].q, .r = 0 };
        }
        res3m[2] = .{ .p = res[4].p, .q = res[4].q, .r = 0 };
    }

    if (comptime family == .x86_64) {
        // AMD parts with AVX-512 run the Skylake-X kernels on other blocking
        const eq = std.mem.eql;
        const skx = comptime (eq(u8, core_name, "skylakex") or eq(u8, core_name, "cooperlake") or eq(u8, core_name, "sapphirerapids"));
        if (skx and env.amd_avx512 and env.l3_kb > 0 and @rem(env.l3_kb, 32768) == 0 and env.l2_kb == 1024) {
            const pq = [4][2]i32{ .{ 384, 512 }, .{ 512, 512 }, .{ 160, 480 }, .{ 176, 256 } };
            inline for (0..4) |i| res[i] = .{ .p = pq[i][0], .q = pq[i][1], .r = 0 };
        }
    }

    // only x86 rounds `p` up to the unroll (`setparam-ref.c:2150-2175`): `s`,
    // `d`, `c`, `z`, and 3M to its own unroll, else the real one's; `q`, `x`
    // and `x3m` only under `QUAD_PRECISION` (`:2179-2183`), off by default
    // (`Makefile.rule:214`)
    if (x86) {
        inline for (.{ "s", "d", "c", "z" }, 0..) |name, i| res[i].p = roundUp(res[i].p, @field(table, name).unroll_m);
        const unroll_of = [2]i32{ table.s.unroll_m, table.d.unroll_m };
        inline for (.{ "c3m", "z3m" }, 0..) |name, i| {
            const um = @field(table, name).unroll_m;
            res3m[i].p = roundUp(res3m[i].p, if (um > 0) um else unroll_of[i]);
        }
    }

    // x86 always derives `r`, the other families take the header's constant
    inline for (names, 0..) |name, i| {
        const fixed = if (x86) null else fixedOf(@field(table, name).r);
        res[i].r = fixed orelse derivedR(res[i].p, res[i].q, elems[i], table);
    }
    if (x86) {
        const elems3m = [3]i64{ 8, 16, 32 };
        inline for (0..3) |i| res3m[i].r = derivedR(res3m[i].p, res3m[i].q, elems3m[i], table);
    }

    // `setparam-ref.c:1099` (`ARCH_MIPS64`), over the header's 8192 or the
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
    inline for (names, 0..) |name, i| setPrec(vt, name, @field(table, name), res[i]);
    inline for (names3m, 0..) |name, i| setPrec(vt, name, @field(table, name), res3m[i]);

    // half and bfloat16: the headers' `p` and `q`; x86 derives `r`
    // (`setparam-ref.c:1556-1563, 2190-2205`), the others take the header's.
    // arm64, LoongArch, POWER and zarch set only `b` and `sb` (`:981-1015,
    // 1181-1186, 1307, 1316-1335, 1348-1367`), riscv64 also `sh`
    // (`:1381-1410`), mips64 none; the static families' are the header's
    // (`common_param.h:1441-1465`). Nothing reads these until their level 3 is
    // ported.
    inline for (.{ "sh", "b", "sb" }) |name| {
        const set = comptime switch (family) {
            .mips64 => false,
            .arm64, .loongarch64, .power64, .power64be, .power32, .power32le, .s390x => !std.mem.eql(u8, name, "sh"),
            .x86_64, .x86, .riscv64, .arm, .mips32, .sparc, .wasm, .generic => true,
        };
        const sh = @field(table, name);
        const p = fixedOf(sh.p) orelse 0;
        const q = fixedOf(sh.q) orelse 0;
        const fixed = if (x86) null else fixedOf(sh.r);
        setPrec(vt, name, sh, if (set) .{ .p = p, .q = q, .r = fixed orelse derivedR(p, q, 4, table) } else .{ .p = 0, .q = 0, .r = 0 });
    }
}

// ---- lookup ----------------------------------------------------------------

/// The architecture families the data is extracted for. The 32-bit ones take
/// `Makefile.system`'s 32-bit fallbacks (x86's AVX cores Nehalem's entry, the
/// Bulldozer family and Zen Barcelona's, POWER8 to 10 POWER6's); POWER's
/// follow the byte order too (`PPC970`'s unrolls).
pub const Family = enum { x86_64, x86, arm64, arm, power64, power64be, power32, power32le, riscv64, s390x, mips64, mips32, loongarch64, sparc, wasm, generic };

/// The core enum of `family`.
pub fn CoreOf(comptime family: Family) type {
    return switch (family) {
        .x86_64, .x86 => coretypes.X86CoreType,
        .arm64 => coretypes.Arm64CoreType,
        .arm => coretypes.ArmCoreType,
        .power64, .power64be, .power32, .power32le => coretypes.PowerCoreType,
        .riscv64 => coretypes.RiscvCoreType,
        .s390x => coretypes.S390xCoreType,
        .mips64, .mips32 => coretypes.MipsCoreType,
        .loongarch64 => coretypes.LoongarchCoreType,
        .sparc => coretypes.SparcCoreType,
        .wasm => coretypes.WasmCoreType,
        .generic => coretypes.GenericCoreType,
    };
}

/// The extracted table of `core`; null where `param.h` has none for it.
pub fn tableFor(comptime family: Family, core: CoreOf(family)) ?Table {
    return switch (family) {
        inline else => |f| @field(data, @tagName(f))(core),
    };
}

/// `DTB_DEFAULT_ENTRIES` of `core` (`getarch.c`). Where `getarch.c` has no
/// entry, a build for the core takes what its `cpuid_<arch>.c` prints (riscv64:
/// 64, for every core), else the generic target's.
pub fn dtbFor(comptime family: Family, core: CoreOf(family)) i32 {
    const own = switch (family) {
        inline else => |f| @field(data, @tagName(f) ++ "_dtb")(core),
    };
    return own orelse switch (family) {
        .riscv64 => 64,
        else => data.generic_dtb(.generic).?,
    };
}

/// Installs the parameters of `core`. A core `param.h` has no blocking for takes
/// the generic set and its rule, as OpenBLAS's generic target does, with its
/// family's `BUFFER_SIZE`.
pub fn installCore(vt: *Vtable, comptime family: Family, core: CoreOf(family), env: Env) void {
    switch (core) {
        inline else => |c| {
            defer vt.dtb_entries = dtbFor(family, c);
            if (tableFor(family, c)) |table| {
                install(vt, family, @tagName(c), table, env);
            } else {
                var table = data.generic(.generic).?;
                table.buffer_size = @field(data, @tagName(family) ++ "_buffer_size");
                install(vt, .generic, "generic", table, env);
            }
        },
    }
}

// ---- tests -----------------------------------------------------------------

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
    inline for (.{ "s", "d", "c", "z" }) |pre| {
        const p = @field(vt, pre ++ "gemm_p");
        const q = @field(vt, pre ++ "gemm_q");
        const r = @field(vt, pre ++ "gemm_r");
        const um = @field(vt, pre ++ "gemm_unroll_m");
        // a precision the core has no header set for stays zero; only x86
        // rounds `p` to the unroll
        const x86 = family == .x86 or family == .x86_64;
        if (um != 0 and (p <= 0 or q <= 0 or r <= 0 or (x86 and @rem(p, um) != 0))) {
            std.debug.print("{s}/{s}/{s}: p={d} q={d} r={d} unroll_m={d}\n", .{ @tagName(family), @tagName(core), pre, p, q, r, um });
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

const q_fields = .{ "sgemm_q", "dgemm_q", "cgemm_q", "zgemm_q", "qgemm_q", "xgemm_q", "cgemm3m_q", "zgemm3m_q", "xgemm3m_q", "shgemm_q", "bgemm_q", "sbgemm_q" };

/// Every precision's `q`.
fn qOf(vt: *const Vtable) [q_fields.len]i32 {
    var out: [q_fields.len]i32 = undefined;
    inline for (q_fields, 0..) |f, i| out[i] = @field(vt, f);
    return out;
}

/// Every core of `family` resolves to sane parameters; one test per family,
/// so each fits a rig image (`build/qemutest.zig`'s `sparc64-debug`).
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

test "generic follows the GENERIC rule for a 512 KB L2" {
    var vt: Vtable = .{};
    installCore(&vt, .generic, .generic, .{});
    // min(512 >> 6, 16) = 8
    try std.testing.expectEqual(@as(i32, 448), vt.sgemm_p);
    try std.testing.expectEqual(@as(i32, 224), vt.dgemm_p);
    try std.testing.expectEqual(@as(i32, 112), vt.zgemm_p);
    try std.testing.expectEqual(@as(i32, 128), vt.dgemm_q);
    const packed_a: i64 = (448 * 128 * 4 + 0xffff) & ~@as(i64, 0xffff);
    const want: i64 = (@divTrunc(data.generic_buffer_size - packed_a, 128 * 4) - 15) & ~@as(i64, 15);
    try std.testing.expectEqual(@as(i32, @intCast(want)), vt.sgemm_r);
}

test "the block factor scales a computed p in 8s, and the L2 size moves it" {
    var vt: Vtable = .{};
    installCore(&vt, .generic, .generic, .{ .block_factor = 50 });
    try std.testing.expectEqual(@as(i32, 224), vt.sgemm_p);
    installCore(&vt, .generic, .generic, .{ .block_factor = 1 });
    // clamped to 10 percent: 44 in 8s is 40
    try std.testing.expectEqual(@as(i32, 40), vt.sgemm_p);
    installCore(&vt, .generic, .generic, .{ .l2_kb = 128 });
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
    // no 3M macros: `setparam-ref.c:1582, 1590, 2097, 2105`
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
    // header's 112 and 224 (`:1597, 2110`)
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
    // the header's 44 on an unroll of 8, which OpenBLAS keeps (`setparam-ref.c:1089`;
    // static mips32 builds take the header's)
    inline for (.{ .mips64, .mips32 }) |f| {
        installCore(&vt, f, .loongson3r4, .{});
        try std.testing.expectEqual(@as(i32, 8), vt.dgemm_unroll_m);
        try std.testing.expectEqual(@as(i32, 44), vt.dgemm_p);
        try std.testing.expectEqual(@as(usize, 40), driver.Blocking.of(f64, &vt).p);
    }
    // core2's rule: 92 * (512 >> 9) + 8 = 100, rounded to 104; `q`'s
    // 23 + 4 = 27 stays (`:2179-2183`)
    installCore(&vt, .x86_64, .core2, .{ .l2_kb = 512 });
    try std.testing.expectEqual(@as(i32, 104), vt.sgemm_p);
    try std.testing.expectEqual(@as(i32, 104), vt.cgemm3m_p);
    try std.testing.expectEqual(@as(i32, 27), vt.qgemm_p);
}

test "half and bfloat16 take the header's r outside x86, and only what each family's block sets" {
    var vt: Vtable = .{};
    // `b` and `sb` only: `setparam-ref.c:1015, 1186, 1325, 1357`
    inline for (.{ .{ .arm64, .neoversen1 }, .{ .loongarch64, .la464 }, .{ .power64, .power10 }, .{ .s390x, .z14 } }) |c| {
        installCore(&vt, c[0], c[1], .{});
        const t = tableFor(c[0], c[1]).?;
        try std.testing.expectEqual(t.sb.r.fixed, vt.sbgemm_r);
        try std.testing.expectEqual(t.b.r.fixed, vt.bgemm_r);
        try std.testing.expectEqual(t.sb.p.fixed, vt.sbgemm_p);
        inline for (.{ "_p", "_q", "_r" }) |f| try std.testing.expectEqual(@as(i32, 0), @field(vt, "shgemm" ++ f));
    }
    installCore(&vt, .arm64, .neoversen1, .{});
    try std.testing.expectEqual(@as(i32, 256), vt.sbgemm_r);
    // riscv64 also `sh`, whose header `r` is 240 (`:1397`)
    installCore(&vt, .riscv64, .riscv64_generic, .{});
    try std.testing.expectEqual(@as(i32, 240), vt.shgemm_r);
    try std.testing.expectEqual(@as(i32, 256), vt.sbgemm_r);
    try std.testing.expectEqual(@as(i32, 256), vt.bgemm_r);
    // mips64 none
    installCore(&vt, .mips64, .loongson3r4, .{});
    inline for (.{ "sh", "b", "sb" }) |pre| {
        inline for (.{ "gemm_p", "gemm_q", "gemm_r" }) |f| try std.testing.expectEqual(@as(i32, 0), @field(vt, pre ++ f));
    }
    // x86 derives it (`:2190-2205`)
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
