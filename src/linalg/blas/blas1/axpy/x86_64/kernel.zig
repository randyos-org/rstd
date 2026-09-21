//! x86_64 axpy kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`).
//!
//! Each core's kernels are the reference's `{s,d,c,z}axpy.c` as compiled for
//! that core: the shared driver (`kernel.zig`) around a block kernel that is
//! either a hand-written inline-asm microkernel, transcribed by script into
//! `microk.zig` and used here as written, or, for the AVX-512 cores
//! whose upstream block is C intrinsics, the same loops on `@Vector`s. Which
//! block a core uses is the reference's `#if defined(<CORE>)` chain.
//!
//! The build gives the object the core's CPU features; the assembler needs
//! them to accept the mnemonics, while whether the C-derived multiply-adds are
//! contracted is the separate policy in `fma.zig`.
//!
//! The cores without a `.c` kernel (Prescott, Core2, Penryn, Atom, Opteron, ...)
//! use whole-function `.S` kernels upstream. Those are in `asm/`, one
//! Zig file per source with its `#if`s and macros resolved for a core and ABI,
//! and are emitted here as naked functions whose body is that assembly; they
//! read their arguments straight from the ABI registers, as the reference's
//! do. The x87 `qaxpy`/`xaxpy` are not ported: they only exist upstream in
//! `EXPRECISION` builds.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../../types.zig");
const coretypes = @import("../../../cpu/_coretypes.zig");
const fma = @import("../../../fma.zig");
const kernel = @import("../kernel.zig");
const mk = @import("microk.zig");
const asmexport = @import("../../../asm/exports.zig");
const axpy_sse = @import("asm/axpy_sse.zig");
const axpy_sse2 = @import("asm/axpy_sse2.zig");
const axpy_atom = @import("asm/axpy_atom.zig");
const zaxpy_sse = @import("asm/zaxpy_sse.zig");
const zaxpy_sse2 = @import("asm/zaxpy_sse2.zig");
const zaxpy_atom = @import("asm/zaxpy_atom.zig");

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

// ---- the whole-function `.S` kernels ---------------------------------------

/// The ABI family whose variant of the `.S` kernels this build uses.
const windows = sysinfo.os.tag == .windows or sysinfo.os.tag == .uefi;

/// `PROLOGUE`'s `.align`, in bytes: 512 for ELF (SysV), 16 for Win64.
const entry_alignment: u29 = if (windows) 16 else 512;

/// The `.S` source `role` (`saxpy`, `daxpy`, `caxpy`, `caxpyc`, `zaxpy`,
/// `zaxpyc`) uses on `core`, if the reference uses one there: its
/// `KERNEL.<CORE>` choice, else the default in `KERNEL`.
fn sFile(comptime role: []const u8, comptime core: Core) ?type {
    const single = std.mem.eql(u8, role, "saxpy");
    const double = std.mem.eql(u8, role, "daxpy");
    const complex_single = std.mem.eql(u8, role, "caxpy") or std.mem.eql(u8, role, "caxpyc");
    return switch (core) {
        .prescott, .core2, .penryn, .dunnington, .opteron, .opteron_sse3, .barcelona, .bobcat, .nano => if (single)
            axpy_sse
        else if (double)
            axpy_sse2
        else if (complex_single)
            zaxpy_sse
        else
            zaxpy_sse2,
        .atom => if (single)
            axpy_sse
        else if (double)
            axpy_atom
        else if (complex_single)
            zaxpy_sse
        else
            zaxpy_atom,
        .nehalem => if (single or double) null else if (complex_single) zaxpy_sse else zaxpy_sse2,
        .bulldozer => if (single) axpy_sse else null,
        else => null,
    };
}

const SText = struct { alignment: u29, code: []const u8 };

/// `role`'s `.S` kernel for `core` on this ABI. The local labels carry the
/// role and core, so the kernels of one object do not clash.
fn sText(comptime role: []const u8, comptime core: Core) SText {
    comptime {
        const conj = role[role.len - 1] == 'c';
        const double = role[0] == 'd' or role[0] == 'z';
        return .{
            .alignment = entry_alignment,
            .code = sFile(role, core).?.text(double, conj, windows, core, asmexport.labelSuffix(role, @tagName(core))),
        };
    }
}

// ---- kernels and exports ---------------------------------------------------

/// Whether `core` has a `role` kernel (`daxpy`, `saxpy`, `caxpy`, `zaxpy`).
pub fn has(comptime core: Core, comptime role: []const u8) bool {
    return comptime blk: {
        if (std.mem.eql(u8, role, "daxpy")) break :blk realBlock(f64, core).block != .none or sFile("daxpy", core) != null;
        if (std.mem.eql(u8, role, "saxpy")) break :blk realBlock(f32, core).block != .none or sFile("saxpy", core) != null;
        if (std.mem.eql(u8, role, "caxpy")) break :blk cplxKernel(f32, core) != null or sFile("caxpy", core) != null;
        break :blk cplxKernel(f64, core) != null or sFile("zaxpy", core) != null;
    };
}

/// Exported symbol of `role` (`daxpy`, `caxpyc`, ...) for `core`, following
/// the reference's `TSUFFIX = _<CORE>` naming.
pub fn symbol(comptime role: []const u8, comptime core: Core) []const u8 {
    return asmexport.symbol(role, @tagName(core));
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

fn exportS(comptime role: []const u8, comptime core: Core) void {
    if (comptime sFile(role, core) != null) @export(&asmexport.Naked(sText(role, core).code, sText(role, core).alignment).f, .{ .name = symbol(role, core) });
}

/// Exports every axpy kernel `core` has. Call from a container-level `comptime`
/// block of the object's root file.
pub fn exportCore(comptime core: Core) void {
    if (has(core, "daxpy")) {
        if (realBlock(f64, core).block != .none) @export(realKernel(f64, core), .{ .name = symbol("daxpy", core) }) else exportS("daxpy", core);
    }
    if (has(core, "saxpy")) {
        if (realBlock(f32, core).block != .none) @export(realKernel(f32, core), .{ .name = symbol("saxpy", core) }) else exportS("saxpy", core);
    }
    if (has(core, "caxpy")) {
        if (cplxKernel(f32, core) != null) {
            @export(cplxEntry(f32, core, false), .{ .name = symbol("caxpy", core) });
            @export(cplxEntry(f32, core, true), .{ .name = symbol("caxpyc", core) });
        } else {
            exportS("caxpy", core);
            exportS("caxpyc", core);
        }
    }
    if (has(core, "zaxpy")) {
        if (cplxKernel(f64, core) != null) {
            @export(cplxEntry(f64, core, false), .{ .name = symbol("zaxpy", core) });
            @export(cplxEntry(f64, core, true), .{ .name = symbol("zaxpyc", core) });
        } else {
            exportS("zaxpy", core);
            exportS("zaxpyc", core);
        }
    }
}
