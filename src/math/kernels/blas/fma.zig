//! Whether a multiply-add is fused, as the reference's per-core compile flags
//! decide it; the policy is the root's `std_float_options` (design.md,
//! Contraction). `linalg.vector.mulAdd` delegates to `mulAdd`.

const std = @import("std");
const sysinfo = @import("builtin");
const root = @import("root");
const types = @import("types.zig");
const util = @import("util.zig");

/// A root's `std_float_options` type, re-exported as `zigcfg.FloatOptions`:
/// `std.Options.Float` where std has it (this project's fork), else a copy
/// std ignores and rstd reads. Here, not in `zigcfg`, so the per-core
/// objects, whose module is this directory, can name it.
pub const FloatOptions = if (has_float_options) std.Options.Float else struct {
    /// Multiply-adds that may fuse do in floats of at most this many bits: 0
    /// never fuses, 128 always does.
    fuse_mul_add_bits: u16 = defaultFuseMulAddBits(),
};

/// Whether the compiler reads `std_float_options` itself, for `@Matrix`.
pub const has_float_options = @hasDecl(std.Options, "Float");

/// The build's float options, the root's or the defaults.
pub const float_options: FloatOptions = if (has_float_options)
    std.float_options
else if (@hasDecl(root, "std_float_options")) root.std_float_options else .{};

const fuse_bits = float_options.fuse_mul_add_bits;
const hardware_bits = fusedMulAddBits(&sysinfo.target);

/// The widest float, in bits, `target` has a fused multiply-add instruction
/// for; 0 for none, or a soft-float ABI. Copied from the fork's
/// `std.Target.fusedMulAddBits` (zig-build/zig-src lib/std/Target.zig,
/// b5896d765b) for stock Zig's `zigcfg.FloatOptions`; re-derive it from there
/// rather than hand-edit.
pub fn fusedMulAddBits(target: *const std.Target) u16 {
    if (target.abi.float() == .soft) return 0;
    const cpu = target.cpu;
    return switch (cpu.arch) {
        .x86, .x86_64 => if (cpu.has(.x86, .fma)) 64 else 0,
        .aarch64, .aarch64_be => 64,
        .arm, .armeb, .thumb, .thumbeb => if (cpu.hasAny(.arm, &.{ .vfp4, .vfp4d16 }))
            64
        else if (cpu.hasAny(.arm, &.{ .vfp4sp, .vfp4d16sp })) 32 else 0,
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => 64,
        .riscv32, .riscv32be, .riscv64, .riscv64be => if (cpu.has(.riscv, .d))
            64
        else if (cpu.has(.riscv, .f)) 32 else 0,
        .loongarch32, .loongarch64 => if (cpu.has(.loongarch, .d))
            64
        else if (cpu.has(.loongarch, .f)) 32 else 0,
        .s390x => 64,
        // before release 6, `madd.fmt` rounds the product
        .mips, .mipsel, .mips64, .mips64el => if (cpu.hasAny(.mips, &.{ .mips32r6, .mips64r6 })) 64 else 0,
        else => 0,
    };
}

/// The fork's default `fuse_mul_add_bits`, for stock Zig's
/// `zigcfg.FloatOptions`: none where LLVM doesn't compile `@mulAdd` right.
pub fn defaultFuseMulAddBits() u16 {
    if (!compiles and sysinfo.zig_backend == .stage2_llvm) return 0;
    return fusedMulAddBits(&sysinfo.target);
}

/// Whether the build fuses at all, for kernels taking one flag for every
/// type; `madd` narrows it per type.
pub const fuse = fuse_bits > 0;

/// Whether `madd.fmt` rounds once, as on Loongson's GS464 and GS464E (their
/// manuals; gcc's `ISA_HAS_FUSED_MADD4`): MIPS before release 6, where the
/// root's `blas_fused_madd` says so. Only the `loongson3r3` and `loongson3r4`
/// objects' does (`cpu/corespec.zig`). LLVM also picks `madd.fmt` for a
/// plain `a*b + c` there, which then fuses too; QEMU rounds it (issues #196, #224).
pub const fused_madd = sysinfo.target.cpu.arch.isMips() and hardware_bits == 0 and
    @hasDecl(root, "blas_fused_madd") and root.blas_fused_madd;

/// Whether `T`, a float or a vector of them, has a fused instruction: up to
/// the hardware's width where `compiles`, or `fused_madd`; never past `f64`.
pub fn fusible(comptime T: type) bool {
    const E = Elem(T);
    return !types.isExtended(E) and
        ((@typeInfo(E).float.bits <= hardware_bits and compiles) or fused_madd);
}

fn Elem(comptime T: type) type {
    return switch (@typeInfo(T)) {
        .vector => |v| v.child,
        else => T,
    };
}

/// Whether the setting lets multiply-adds in `T` fuse.
fn allowed(comptime T: type) bool {
    return @typeInfo(Elem(T)).float.bits <= fuse_bits;
}

/// Whether LLVM compiles `@mulAdd` right: not with MSA, where a vector `fma`,
/// also one the vectorizers make of scalar ones, becomes `fmadd` with the
/// addend as a multiplicand (issue #194).
pub const compiles = !(sysinfo.target.cpu.arch.isMips() and std.Target.mips.featureSetHas(sysinfo.target.cpu.features, .msa));

/// Whether the BLAS fuses multiply-adds in `T`: the setting allows it and an
/// instruction does it.
pub fn fuseFor(comptime T: type) bool {
    return allowed(T) and fusible(T);
}

/// `a*b + c`, fused if `fused` and `fuseFor`; kept out of AltiVec's flushing
/// lanes (`util.scalar`). For kernels that pick their own fusing (`fuse`).
pub inline fn madd(comptime fused: bool, comptime T: type, a: T, b: T, c: T) T {
    return util.scalar(T, if (comptime fused and fuseFor(T)) fusedMadd(T, a, b, c) else a * b + c);
}

/// `fuseFor` and `allowed` by float width, for `mulAdd`: worked out here
/// once, since a comptime call charges its whole branch count to every
/// inline caller, memoized or not.
const by_bits = blk: {
    var r: [129]struct { fuse: bool, allow: bool } = @splat(.{ .fuse = false, .allow = false });
    for (.{ f16, f32, f64, f80, f128 }) |F| r[@typeInfo(F).float.bits] = .{ .fuse = fuseFor(F), .allow = allowed(F) };
    break :blk r;
};

/// `a*b + c` by the setting alone, as the fork's `@Matrix` products: fused
/// in floats it allows, in software past the hardware's width.
pub inline fn mulAdd(comptime T: type, a: T, b: T, c: T) T {
    const E = switch (@typeInfo(T)) {
        .vector => |v| v.child,
        else => T,
    };
    const policy = comptime by_bits[@typeInfo(E).float.bits];
    if (policy.fuse) return fusedMadd(T, a, b, c);
    return if (policy.allow and compiles) @mulAdd(T, a, b, c) else a * b + c;
}

/// `a*b + c` rounded once, for a `fusible` `T`: `@mulAdd`, or with
/// `fused_madd` `madd.fmt` per scalar and MSA's `fmadd` per 128 bits (its
/// `wd + ws*wt`; LLVM's own `fmadd` swaps operands, issue #194).
pub inline fn fusedMadd(comptime T: type, a: T, b: T, c: T) T {
    if (comptime !fused_madd) return @mulAdd(T, a, b, c);
    const v = switch (@typeInfo(T)) {
        .vector => |v| v,
        else => return asm ((if (T == f32) "madd.s" else "madd.d") ++ " %[d], %[c], %[a], %[b]"
            : [d] "=f" (-> T),
            : [c] "f" (c),
              [a] "f" (a),
              [b] "f" (b),
        ),
    };
    const E = v.child;
    const w = 16 / @sizeOf(E);
    const msa = comptime std.Target.mips.featureSetHas(sysinfo.target.cpu.features, .msa);
    if (msa and v.len == w) return asm ((if (E == f32) "fmadd.w" else "fmadd.d") ++ " %[d], %[a], %[b]"
        : [d] "=f" (-> T),
        : [c] "0" (c),
          [a] "f" (a),
          [b] "f" (b),
    );
    var r: T = undefined;
    if (msa and v.len % w == 0) {
        // whole registers, taken out and put back lane by lane
        inline for (0..v.len / w) |k| {
            const lanes = comptime std.simd.iota(i32, w) + @as(@Vector(w, i32), @splat(k * w));
            const W = @Vector(w, E);
            const part = fusedMadd(W, @shuffle(E, a, undefined, lanes), @shuffle(E, b, undefined, lanes), @shuffle(E, c, undefined, lanes));
            inline for (0..w) |j| r[k * w + j] = part[j];
        }
    } else inline for (0..v.len) |i| {
        r[i] = fusedMadd(E, a[i], b[i], c[i]);
    }
    return r;
}

test "mulAdd agrees with a*b+c on exact values" {
    try std.testing.expectEqual(@as(f64, 14.0), mulAdd(f64, 2.0, 3.0, 8.0));
    try std.testing.expectEqual(@as(f32, 14.0), mulAdd(f32, 2.0, 3.0, 8.0));
}

test "fusedMulAddBits matches the fork's" {
    if (@hasDecl(std.Target, "fusedMulAddBits")) {
        try std.testing.expectEqual(sysinfo.target.fusedMulAddBits(), fusedMulAddBits(&sysinfo.target));
        try std.testing.expectEqual((std.Options.Float{}).fuse_mul_add_bits, defaultFuseMulAddBits());
    } else return error.SkipZigTest;
}

test "extended floats never fuse" {
    inline for (.{ f80, f128, @Vector(2, f80) }) |T| try std.testing.expect(!fuseFor(T));
    // by width, not identity: a 64-bit `c_longdouble` of its own fuses as f64 does
    try std.testing.expectEqual(fuseFor(f64) and @typeInfo(c_longdouble).float.bits <= 64, fuseFor(c_longdouble));
    // (1 + 2^-63)(1 - 2^-63) - 1 is -2^-126 fused, 0 rounded first
    const e: f80 = std.math.ldexp(@as(f80, 1), -63);
    try std.testing.expectEqual(@as(f80, 0), mulAdd(f80, 1 + e, 1 - e, -1));
}
