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
//!   * `p` of the older x86 cores follows the L2 size (`Affine`);
//!   * on x86, `r` is always derived from the packing buffer;
//!   * `p` is rounded up to the unroll;
//!   * the 3M variants take `p` and `q` of the real ones unless the header says;
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

/// How `p`, `q` or `r` comes about.
pub const Size = union(enum) {
    /// A constant of `param.h`.
    fixed: i32,
    /// The header names a run-time variable; see `Rule`.
    dynamic,
    /// `r` derived from the packing buffer once `p` and `q` are known.
    from_buffer,
};

/// One precision's parameters.
pub const Shape = struct {
    unroll_m: i32 = 0,
    unroll_n: i32 = 0,
    /// `GEMM_UNROLL_MN`.
    unroll_mn: i32 = 0,
    /// `BGEMM_ALIGN_K` and `SBGEMM_ALIGN_K`.
    align_k: i32 = 0,
    /// `<P>GEMM_DEFAULT_P_BASE`, where the header has one.
    p_base: i32 = 0,
    p: Size = .{ .fixed = 0 },
    q: Size = .{ .fixed = 0 },
    r: Size = .{ .fixed = 0 },
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
};

/// The packing buffer `r` is derived from (`BUFFER_SIZE`, `common_<arch>.h`).
pub const buffer_size: i64 = 32 << 20;

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
    switch (sysinfo.cpu.arch) {
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
    switch (sysinfo.cpu.arch) {
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
    switch (sysinfo.cpu.arch) {
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
    if (comptime family == .loongarch64 and eq(u8, core_name, "la464")) return .la464;
    if (comptime family == .loongarch64 and eq(u8, core_name, "la264")) return .la264;
    if (comptime family != .x86 and family != .x86_64) return .none;
    if (comptime eq(u8, core_name, "katmai") or eq(u8, core_name, "coppermine") or eq(u8, core_name, "banias"))
        return .{ .affine = .{ .mul = .{ 64, 32, 32, 16, 16, 8 } } };
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
    const v: i64 = @intFromFloat(@as(f64, @floatFromInt(p)) * factor * 1.0e-2);
    return @intCast(v & ~@as(i64, 7));
}

fn fixedOf(s: Size) ?i32 {
    return switch (s) {
        .fixed => |v| v,
        else => null,
    };
}

/// `p` of one of the primary precisions (`idx`: s, d, c, z, q, x).
fn primaryP(shape: Shape, rule: Rule, idx: usize, env: Env) i32 {
    const p: i32 = if (fixedOf(shape.p)) |v| v else switch (rule) {
        .affine => |a| blk: {
            const size = if (a.cap > 0) @min(env.l2_kb >> a.shift, a.cap) else env.l2_kb >> a.shift;
            break :blk scaled(a.mul[idx] * size + a.add[idx], env);
        },
        else => 0,
    };
    return roundUp(p, shape.unroll_m);
}

/// `r` from the buffer: `setparam-ref.c`'s formula.
fn derivedR(p: i32, q: i32, elem: i64, table: Table) i32 {
    if (q <= 0) return 0;
    const align_: i64 = table.align_;
    const packed_a = (@as(i64, p) * q * elem + table.offset_a + align_) & ~align_;
    const r = @divTrunc(buffer_size - packed_a, @as(i64, q) * elem) - 15;
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

    // the primaries, whose `p` the rules and overrides decide; x86 always
    // derives `r`, the other families take the header's constant
    const names = .{ "s", "d", "c", "z", "q", "x" };
    const elems = [6]i64{ 4, 8, 8, 16, 16, 32 };
    var res: [6]Resolved = undefined;
    inline for (names, 0..) |name, i| {
        const sh = @field(table, name);
        const p = primaryP(sh, rule, i, env);
        const q = fixedOf(sh.q) orelse 0;
        const r = if (x86) derivedR(p, q, elems[i], table) else (fixedOf(sh.r) orelse derivedR(p, q, elems[i], table));
        res[i] = .{ .p = p, .q = q, .r = r };
    }

    if (comptime std.meta.activeTag(rule) == .la464) {
        const set = la464_sets[@intFromBool(env.threads == 1)][@intFromBool(@divTrunc(env.l3_kb, 1024) == 32)];
        inline for (.{ "s", "d", "c", "z" }, 0..) |name, i| {
            res[i] = .{ .p = roundUp(set[i].p, @field(table, name).unroll_m), .q = set[i].q, .r = set[i].r };
        }
    }
    if (comptime std.meta.activeTag(rule) == .la264) {
        // the 2K3000 (processor id 0x0014b020) has its own `z` set
        const set: Resolved = if (env.prid == 0x0014b020) .{ .p = 128, .q = 176, .r = 360 } else .{ .p = 64, .q = 120, .r = 4096 };
        res[3] = .{ .p = roundUp(set.p, table.z.unroll_m), .q = set.q, .r = set.r };
    }
    if (comptime x86) {
        // AMD parts with AVX-512 run the Skylake-X kernels on other blocking
        const eq = std.mem.eql;
        const skx = comptime (eq(u8, core_name, "skylakex") or eq(u8, core_name, "cooperlake") or eq(u8, core_name, "sapphirerapids"));
        if (skx and env.amd_avx512 and env.l3_kb > 0 and @rem(env.l3_kb, 32768) == 0 and env.l2_kb == 1024) {
            const pq = [4][2]i32{ .{ 384, 512 }, .{ 512, 512 }, .{ 160, 480 }, .{ 176, 256 } };
            inline for (0..4) |i| {
                res[i].p = pq[i][0];
                res[i].q = pq[i][1];
                res[i].r = derivedR(res[i].p, res[i].q, elems[i], table);
            }
        }
    }
    inline for (names, 0..) |name, i| setPrec(vt, name, @field(table, name), res[i]);

    // half and bfloat16: their headers' `p` and `q`, `r` from the buffer
    inline for (.{ "sh", "b", "sb" }) |name| {
        const sh = @field(table, name);
        const p = roundUp(fixedOf(sh.p) orelse 0, sh.unroll_m);
        const q = fixedOf(sh.q) orelse 0;
        setPrec(vt, name, sh, .{ .p = p, .q = q, .r = derivedR(p, q, 4, table) });
    }

    // 3M: the header's `p` and `q`, else those of the real precision
    const from = [3]usize{ 0, 1, 4 };
    const elem_of = [3]usize{ 2, 3, 5 };
    const unroll_of = [3]i32{ table.s.unroll_m, table.d.unroll_m, table.q.unroll_m };
    inline for (.{ "c3m", "z3m", "x3m" }, 0..) |name, i| {
        const sh = @field(table, name);
        const p = roundUp(fixedOf(sh.p) orelse res[from[i]].p, if (sh.unroll_m > 0) sh.unroll_m else unroll_of[i]);
        const q = fixedOf(sh.q) orelse res[from[i]].q;
        setPrec(vt, name, sh, .{ .p = p, .q = q, .r = derivedR(p, q, elems[elem_of[i]], table) });
    }
}

// ---- lookup ----------------------------------------------------------------

/// The architecture families the data is extracted for.
pub const Family = enum { x86_64, x86, arm64, arm, power64, power32, riscv64, s390x, mips64, mips32, loongarch64, sparc, wasm, generic };

/// The core enum of `family`.
pub fn CoreOf(comptime family: Family) type {
    return switch (family) {
        .x86_64, .x86 => coretypes.X86CoreType,
        .arm64 => coretypes.Arm64CoreType,
        .arm => coretypes.ArmCoreType,
        .power64, .power32 => coretypes.PowerCoreType,
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
/// the generic set and its rule, as OpenBLAS's generic target does.
pub fn installCore(vt: *Vtable, comptime family: Family, core: CoreOf(family), env: Env) void {
    switch (core) {
        inline else => |c| {
            defer vt.dtb_entries = dtbFor(family, c);
            if (tableFor(family, c)) |table| {
                install(vt, family, @tagName(c), table, env);
            } else {
                install(vt, .generic, "generic", data.generic(.generic).?, env);
            }
        },
    }
}

// ---- tests -----------------------------------------------------------------

test "dtb entries: getarch's per core, cpuid's or generic's otherwise" {
    try std.testing.expectEqual(@as(i32, 256), dtbFor(.x86_64, .core2));
    try std.testing.expectEqual(@as(i32, 40), dtbFor(.x86_64, .bobcat));
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
        // a precision the core has no header set for stays zero
        if (um != 0 and (p <= 0 or q <= 0 or r <= 0 or @rem(p, um) != 0)) {
            std.debug.print("{s}/{s}/{s}: p={d} q={d} r={d} unroll_m={d}\n", .{ @tagName(family), @tagName(core), pre, p, q, r, um });
            return error.TestUnexpectedResult;
        }
    }
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

test "every core of the power32 family resolves to sane parameters" {
    try checkFamily(.power32);
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
    const want: i64 = (@divTrunc(buffer_size - packed_a, 128 * 4) - 15) & ~@as(i64, 15);
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

test "old x86 cores follow the L2 size" {
    var vt: Vtable = .{};
    installCore(&vt, .x86_64, .core2, .{ .l2_kb = 4096 });
    // 92 * (4096 >> 9) + 8 = 744
    try std.testing.expectEqual(@as(i32, 744), vt.sgemm_p);
    installCore(&vt, .x86_64, .atom, .{ .l2_kb = 4096 });
    try std.testing.expectEqual(@as(i32, 256), vt.sgemm_p);
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
