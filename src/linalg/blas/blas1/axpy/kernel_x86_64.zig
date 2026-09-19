//! x86_64 axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`).
//!
//! Each core's kernels are the reference's `{s,d,c,z}axpy.c` as compiled for
//! that core: the shared driver (`kernel.zig`) around a block kernel that is
//! either a hand-written inline-asm microkernel, transcribed by script into
//! `x86_64_microk.zig` and used here as written, or, for the AVX-512 cores
//! whose upstream block is C intrinsics, the same loops on `@Vector`s. Which
//! block a core uses is the reference's `#if defined(<CORE>)` chain.
//!
//! The build gives the object the core's CPU features; the assembler needs
//! them to accept the mnemonics, while whether the C-derived multiply-adds are
//! contracted is the separate policy in `fma.zig`.
//!
//! The cores without a `.c` kernel (Prescott, Core2, Penryn, Atom, Opteron, ...)
//! use whole-function `.S` kernels upstream. Those are preprocessed per core
//! and ABI by script (`x86_64_sfiles.zig`) and emitted here as naked functions
//! whose body is that assembly; they read their arguments straight from the
//! ABI registers, as the reference's do. The x87 `qaxpy`/`xaxpy` are not
//! ported: they only exist upstream in `EXPRECISION` builds.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../types.zig");
const coretypes = @import("../../cpu/_coretypes.zig");
const fma = @import("../../fma.zig");
const kernel = @import("kernel.zig");
const mk = @import("x86_64_microk.zig");
// const sf = @import("x86_64_sfiles.zig"); // .S kernels, disabled for now

const BlasLong = types.BlasLong;
const Complex = std.math.Complex;
pub const Core = coretypes.X86CoreType;

const Clobbers = std.lang.assembly.Clobbers;
const has_avx = std.Target.x86.featureSetHas(sysinfo.cpu.features, .avx);
const has_avx512f = std.Target.x86.featureSetHas(sysinfo.cpu.features, .avx512f);

/// The C clobber list of `K`, widened to the vector registers this build's
/// target actually has (the C lists name `xmmN` while using `ymmN`).
fn clobberSet(comptime K: type) Clobbers {
    comptime {
        var c: Clobbers = .{};
        for (K.clobbers) |name| {
            if (std.mem.eql(u8, name, "cc")) {
                c.cc = true;
            } else if (std.mem.eql(u8, name, "memory")) {
                c.memory = true;
            } else if (std.mem.startsWith(u8, name, "xmm")) {
                const wide = if (has_avx512f) "zmm" else if (has_avx) "ymm" else "xmm";
                @field(c, wide ++ name[3..]) = true;
            } else {
                @field(c, name) = true;
            }
        }
        return c;
    }
}

// ---- block kernels: transcribed asm ----------------------------------------

inline fn blockReal(comptime K: type, comptime text: []const u8, n: usize, x: anytype, y: anytype, a: anytype) void {
    var i: usize = 0;
    var nn = n;
    asm volatile (text
        : [i] "+r" (i),
          [n] "+r" (nn),
        : [x] "r" (x),
          [y] "r" (y),
          [a] "r" (a),
        : clobberSet(K));
}

inline fn blockCplx(comptime K: type, comptime text: []const u8, n: usize, x: anytype, y: anytype, a: anytype, m: anytype) void {
    var i: usize = 0;
    var nn = n;
    asm volatile (text
        : [i] "+r" (i),
          [n] "+r" (nn),
        : [x] "r" (x),
          [y] "r" (y),
          [a] "r" (a),
          [m] "r" (m),
        : clobberSet(K));
}

/// Real block kernel over the transcribed asm `K`; picks between its two
/// blocks at `K.threshold` where the reference's `if (n < N)` does.
fn RealAsm(comptime T: type, comptime K: type) type {
    return struct {
        pub fn run(n1: BlasLong, alpha: T, x: [*]T, y: [*]T) void {
            var a = alpha;
            const n: usize = @intCast(n1);
            if (comptime K.threshold) |th| {
                if (n1 < th) blockReal(K, K.text[0], n, x, y, &a) else blockReal(K, K.text[1], n, x, y, &a);
            } else {
                blockReal(K, K.text[0], n, x, y, &a);
            }
        }
    };
}

fn CplxAsm(comptime F: type, comptime K: type, comptime conj: bool) type {
    const C = Complex(F);
    const texts = if (conj) K.text_conj else K.text;
    const mv = if (conj) K.mvec_conj else K.mvec;
    return struct {
        pub fn run(n1: BlasLong, alpha: C, x: [*]C, y: [*]C) void {
            var a = alpha;
            var m: [mv.len]F = undefined;
            inline for (mv, 0..) |v, k| m[k] = @floatCast(v);
            const n: usize = @intCast(n1);
            const px: [*]F = @ptrCast(x);
            const py: [*]F = @ptrCast(y);
            if (comptime K.threshold) |th| {
                if (n1 < th) blockCplx(K, texts[0], n, px, py, &a, &m) else blockCplx(K, texts[1], n, px, py, &a, &m);
            } else {
                blockCplx(K, texts[0], n, px, py, &a, &m);
            }
        }
    };
}

// ---- block kernels: the AVX-512 intrinsics loops ---------------------------

/// `{s,d}axpy_microk_skylakex-2.c`: `n & ~31` (`~63` for `f32`) elements in
/// 512-bit blocks, then the rest in 256-bit blocks, each four vectors wide.
fn SkylakexBlock(comptime T: type) type {
    const wide = 64 / @sizeOf(T); // lanes in 512 bits
    const narrow = 32 / @sizeOf(T);
    const V5 = @Vector(wide, T);
    const V2 = @Vector(narrow, T);
    return struct {
        pub fn run(n1: BlasLong, alpha: T, x: [*]T, y: [*]T) void {
            const n: usize = @intCast(n1);
            const a5: V5 = @splat(alpha);
            const a2: V2 = @splat(alpha);
            var i: usize = 0;
            const n512 = n & ~@as(usize, 4 * wide - 1);
            while (i < n512) : (i += 4 * wide) {
                inline for (0..4) |k| {
                    const o = i + k * wide;
                    const vx: V5 = x[o..][0..wide].*;
                    const vy: V5 = y[o..][0..wide].*;
                    y[o..][0..wide].* = @mulAdd(V5, a5, vx, vy);
                }
            }
            while (i < n) : (i += 4 * narrow) {
                inline for (0..4) |k| {
                    const o = i + k * narrow;
                    const vx: V2 = x[o..][0..narrow].*;
                    const vy: V2 = y[o..][0..narrow].*;
                    y[o..][0..narrow].* = @mulAdd(V2, a2, vx, vy);
                }
            }
        }
    };
}

// ---- which block a core uses (`{s,d,c,z}axpy.c`'s `#if defined(<CORE>)`) ---

const Block = enum { none, asm_block, skylakex };

fn realBlock(comptime T: type, comptime core: Core) struct { block: Block, K: type } {
    const d = T == f64;
    return switch (core) {
        .nehalem => .{ .block = .asm_block, .K = if (d) mk.daxpy_nehalem else mk.saxpy_nehalem },
        .sandybridge => .{ .block = .asm_block, .K = if (d) mk.daxpy_sandy else mk.saxpy_sandy },
        .haswell, .zen => .{ .block = .asm_block, .K = if (d) mk.daxpy_haswell else mk.saxpy_haswell },
        .bulldozer => if (d) .{ .block = .asm_block, .K = mk.daxpy_bulldozer } else .{ .block = .none, .K = void },
        .piledriver => .{ .block = .asm_block, .K = if (d) mk.daxpy_piledriver else mk.saxpy_piledriver },
        .steamroller, .excavator => .{ .block = .asm_block, .K = if (d) mk.daxpy_steamroller else mk.saxpy_piledriver },
        .skylakex, .cooperlake, .sapphirerapids => .{ .block = .skylakex, .K = void },
        else => .{ .block = .none, .K = void },
    };
}

fn cplxKernel(comptime F: type, comptime core: Core) ?type {
    const d = F == f64;
    return switch (core) {
        .sandybridge => if (d) mk.zaxpy_sandy else mk.caxpy_sandy,
        .haswell, .zen, .skylakex, .cooperlake, .sapphirerapids => if (d) mk.zaxpy_haswell else mk.caxpy_haswell,
        .bulldozer => if (d) mk.zaxpy_bulldozer else mk.caxpy_bulldozer,
        .piledriver, .steamroller, .excavator => if (d) mk.zaxpy_steamroller else mk.caxpy_steamroller,
        else => null,
    };
}

// ---- kernels and exports ---------------------------------------------------

// The whole-function `.S` kernels are disabled for now; see `x86_64_asm.zig`.
// // ---- the whole-function `.S` kernels ---------------------------------------
//
// /// The ABI family whose variant of the `.S` kernels this build uses.
// const abi_os: sf.Os = if (sysinfo.os.tag == .windows or sysinfo.os.tag == .uefi) .win64 else .sysv;
//
// const SText = struct { alignment: u29, code: []const u8 };
//
// /// Joins a text's lines (see `x86_64_sfiles.zig`) into assembly.
// fn join(comptime lines: []const []const u8, comptime idx: []const u16) []const u8 {
//     comptime {
//         @setEvalBranchQuota(10_000_000);
//         var n: usize = 0;
//         for (idx) |i| n += lines[i].len + 1;
//         var buf: [n]u8 = undefined;
//         var p: usize = 0;
//         for (idx) |i| {
//             @memcpy(buf[p..][0..lines[i].len], lines[i]);
//             p += lines[i].len;
//             buf[p] = '\n';
//             p += 1;
//         }
//         const out = buf;
//         return &out;
//     }
// }
//
// /// Index of the `.S` text `role` (`saxpy`, `caxpyc`, ...) uses for `core`, if
// /// the reference uses one there. Cheap: it only looks at the entry table.
// fn sIndex(comptime role: []const u8, comptime core: Core) ?usize {
//     comptime {
//         for (@field(sf, role).entries) |e| {
//             if (e.os == abi_os and std.mem.eql(u8, e.core, @tagName(core))) return e.text;
//         }
//         return null;
//     }
// }
//
// /// Text number `t` of `role`. Memoized per text, so cores that share one
// /// build it once.
// fn sText(comptime role: []const u8, comptime t: usize) SText {
//     comptime {
//         const R = @field(sf, role);
//         return .{ .alignment = R.texts[t].alignment, .code = join(&R.lines, R.texts[t].lines) };
//     }
// }
//
// /// A function whose body is `text` and nothing else; the exported symbol is
// /// referenced through a C-ABI pointer type, so its declared signature is moot.
// fn Naked(comptime t: SText) type {
//     return struct {
//         fn f() align(t.alignment) callconv(.naked) void {
//             asm volatile (t.code);
//         }
//     };
// }

// ---- kernels and exports ---------------------------------------------------

/// Whether `core` has a `role` kernel (`daxpy`, `saxpy`, `caxpy`, `zaxpy`).
pub fn has(comptime core: Core, comptime role: []const u8) bool {
    return comptime blk: {
        if (std.mem.eql(u8, role, "daxpy")) break :blk realBlock(f64, core).block != .none;
        if (std.mem.eql(u8, role, "saxpy")) break :blk realBlock(f32, core).block != .none;
        if (std.mem.eql(u8, role, "caxpy")) break :blk cplxKernel(f32, core) != null;
        break :blk cplxKernel(f64, core) != null;
    };
}

/// Exported symbol of `role` (`daxpy`, `caxpyc`, ...) for `core`, following
/// the reference's `TSUFFIX = _<CORE>` naming.
pub fn symbol(comptime role: []const u8, comptime core: Core) []const u8 {
    comptime {
        var up: [@tagName(core).len]u8 = undefined;
        for (@tagName(core), 0..) |ch, k| up[k] = std.ascii.toUpper(ch);
        return "rstd_blas_" ++ role ++ "_k_" ++ &up;
    }
}

fn realKernel(comptime T: type, comptime core: Core) kernel.RealKernel(T) {
    const r = realBlock(T, core);
    const mask = 128 / @sizeOf(T);
    const B = switch (r.block) {
        .asm_block => RealAsm(T, r.K),
        .skylakex => SkylakexBlock(T),
        .none => null,
    };
    return kernel.realEntry(T, kernel.Real(T, mask, B, fma.fuse));
}

fn cplxEntry(comptime F: type, comptime core: Core, comptime conj: bool) kernel.ComplexKernel(F) {
    const K = cplxKernel(F, core).?;
    return kernel.complexEntry(F, kernel.Cplx(F, conj, 512 / (@sizeOf(F) * 4), CplxAsm(F, K, conj), fma.fuse));
}

// fn exportS(comptime role: []const u8, comptime core: Core) void {
//     if (comptime sIndex(role, core)) |t| @export(&Naked(sText(role, t)).f, .{ .name = symbol(role, core) });
// }

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    if (has(core, "daxpy")) {
        if (realBlock(f64, core).block != .none) @export(realKernel(f64, core), .{ .name = symbol("daxpy", core) });
    }
    if (has(core, "saxpy")) {
        if (realBlock(f32, core).block != .none) @export(realKernel(f32, core), .{ .name = symbol("saxpy", core) });
    }
    if (has(core, "caxpy")) {
        if (cplxKernel(f32, core) != null) {
            @export(cplxEntry(f32, core, false), .{ .name = symbol("caxpy", core) });
            @export(cplxEntry(f32, core, true), .{ .name = symbol("caxpyc", core) });
        }
    }
    if (has(core, "zaxpy")) {
        if (cplxKernel(f64, core) != null) {
            @export(cplxEntry(f64, core, false), .{ .name = symbol("zaxpy", core) });
            @export(cplxEntry(f64, core, true), .{ .name = symbol("zaxpyc", core) });
        }
    }
}
