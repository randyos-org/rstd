const std = @import("std");
const sysinfo = @import("builtin");
const Vtable = @import("Vtable.zig");
const types = @import("types.zig");

const isComplex = types.isComplex;
const Float = types.Float;
const Allocator = std.mem.Allocator;

/// A vtable slot's name; too many fields for the default branch quota.
pub const Slot = blk: {
    @setEvalBranchQuota(10_000);
    break :blk std.meta.FieldEnum(Vtable);
};

/// The function type of a kernel slot.
pub fn SlotFn(comptime slot: Slot) type {
    return @typeInfo(@FieldType(Vtable, @tagName(slot))).optional.child;
}

/// C's `base[offset]` for a negative `offset` too; BLAS strides may be
/// negative, which Zig's `usize` indexing does not take.
pub inline fn at(comptime T: type, base: [*]T, offset: i64) *T {
    const addr: usize = @intFromPtr(base);
    const byte_offset: i64 = offset * @as(i64, @sizeOf(T));
    const new_addr: usize = if (byte_offset >= 0)
        addr + @as(usize, @intCast(byte_offset))
    else
        addr - @as(usize, @intCast(-byte_offset));
    return @ptrFromInt(new_addr);
}

/// `base + offset` for any sign of `offset`; `at` as a many-item pointer.
pub inline fn off(comptime T: type, base: [*]T, offset: i64) [*]T {
    return @ptrCast(at(T, base, offset));
}

/// Where a vector of `n` with stride `inc` starts: its last element when
/// `inc < 0`. The vector must pass `fits`.
pub fn first(comptime T: type, p: [*]T, n: isize, inc: isize) [*]T {
    return if (inc < 0) off(T, p, -(n - 1) * inc) else p;
}

/// `T`'s elements as the kernels take them: floats, two per complex. Drops
/// `const`; the kernels take mutable pointers.
pub fn raw(comptime T: type, p: [*]const T) [*]Float(T) {
    return @ptrCast(@constCast(p));
}

pub fn zero(comptime T: type) T {
    return if (comptime isComplex(T)) .init(0, 0) else 0;
}

pub fn one(comptime T: type) T {
    return fromReal(T, 1);
}

/// `v` as a `T`, with a zero imaginary part.
pub fn fromReal(comptime T: type, v: Float(T)) T {
    return if (comptime isComplex(T)) .init(v, 0) else v;
}

pub fn isZero(comptime T: type, v: T) bool {
    return if (comptime isComplex(T)) v.re == 0 and v.im == 0 else v == 0;
}

pub fn isOne(comptime T: type, v: T) bool {
    return if (comptime isComplex(T)) v.re == 1 and v.im == 0 else v == 1;
}

/// `MAX_STACK_ALLOC`: up to this many bytes of a routine's scratch live on
/// the stack; the rest comes from the allocator.
pub const stack_bytes = 2048;

/// Scratch for `F` elements: in the caller's `stack` when `stack_bytes` holds
/// them, else from the allocator, aligned to `heap_align`.
pub fn Scratch(comptime F: type, comptime stack_align: usize, comptime heap_align: usize) type {
    return struct {
        const Self = @This();
        ptr: [*]F,
        len: usize,
        heap: ?[]align(heap_align) F = null,

        pub fn init(allocator: Allocator, stack: *align(stack_align) [stack_bytes]u8, len: usize) Allocator.Error!Self {
            if (len * @sizeOf(F) <= stack_bytes) {
                // tests poison it as their allocators do the heap: scratch
                // read before it is written shows as NaN
                if (sysinfo.is_test) @memset(stack, 0xff);
                return .{ .ptr = @ptrCast(stack), .len = len };
            }
            return heapOnly(allocator, len);
        }

        /// Always from the allocator, as the reference's `blas_memory_alloc`.
        pub fn heapOnly(allocator: Allocator, len: usize) Allocator.Error!Self {
            const h = try allocator.alignedAlloc(F, .fromByteUnits(heap_align), len);
            return .{ .ptr = h.ptr, .len = len, .heap = h };
        }

        pub fn slice(self: Self) []F {
            return self.ptr[0..self.len];
        }

        pub fn deinit(self: Self, allocator: Allocator) void {
            if (self.heap) |h| allocator.free(h);
        }
    };
}

/// A scalar as the kernels take it, for `@call`: `.{v}`, or `.{ v.re, v.im }`
/// for a complex `T`.
pub inline fn scalarArgs(comptime T: type, v: T) if (isComplex(T)) struct { Float(T), Float(T) } else struct { T } {
    return if (comptime isComplex(T)) .{ v.re, v.im } else .{v};
}

/// Refuses, at compile time, a `T` that is not a float or its `Complex`.
pub fn requireFloat(comptime T: type, comptime name: []const u8) void {
    if (@typeInfo(Float(T)) != .float) @compileError(name ++ ": " ++ @typeName(T) ++ " is not a float type or its Complex");
}

/// Refuses, at compile time, a real `T`.
pub fn requireComplex(comptime T: type, comptime name: []const u8) void {
    if (!isComplex(T)) @compileError(name ++ ": " ++ @typeName(T) ++ " is not complex");
}

/// Refuses, at compile time, a complex `T`.
pub fn requireReal(comptime T: type, comptime name: []const u8) void {
    if (isComplex(T)) @compileError(name ++ ": " ++ @typeName(T) ++ " is not real");
}

/// `@abs(v)` of a float or float vector; on MIPS MSA by a mask LLVM cannot
/// see, as its `@abs` there keeps the sign (`issues.md`). For any `|x|` MSA may vectorize.
pub fn abs(v: anytype) @TypeOf(v) {
    const T = @TypeOf(v);
    const F = if (@typeInfo(T) == .vector) @typeInfo(T).vector.child else T;
    if (comptime !msa or (F != f32 and F != f64)) return @abs(v);
    const lanes = comptime if (@typeInfo(T) == .vector) @typeInfo(T).vector.len else 1;
    const words = @sizeOf(F) / 4;
    const n = lanes * words;
    const k = asm (""
        : [k] "=r" (-> u32),
        : [x] "0" (@as(u32, 0x7fff_ffff)),
    );
    // the sign's word of each float: its last, or on big-endian its first
    const top = comptime if (sysinfo.target.cpu.arch.endian() == .little) words - 1 else 0;
    var mask: [n]u32 = @splat(0xffff_ffff);
    inline for (0..n) |i| {
        if (i % words == top) mask[i] = k;
    }
    return @bitCast(@as(@Vector(n, u32), @bitCast(v)) & @as(@Vector(n, u32), mask));
}

/// Whether this target is x86 with `feature`.
pub fn x86Has(comptime feature: std.Target.x86.Feature) bool {
    return sysinfo.target.cpu.arch.isX86() and std.Target.x86.featureSetHas(sysinfo.target.cpu.features, feature);
}

const msa = sysinfo.target.cpu.arch.isMIPS() and std.Target.mips.featureSetHas(sysinfo.target.cpu.features, .msa);

/// `core`'s kernel symbol for `slot` (`blaskern_daxpy_k_haswell`), exported hidden by
/// `build/blaskern.zig`'s objects and looked up by `kernels/<arch>.zig`'s `install`.
pub fn kernelSymbol(comptime slot: Slot, comptime core: anytype) []const u8 {
    return "blaskern_" ++ @tagName(slot) ++ "_" ++ @tagName(core);
}

/// Loads a `V` from `p`, unaligned; a caller's buffer at any stride need not
/// have `@Vector`'s alignment. Here, not in `simd/`: the per-core objects
/// (`_kernel.zig`) cannot reach outside `linalg/blas/`.
pub inline fn load(comptime V: type, p: [*]const @typeInfo(V).vector.child) V {
    return @as(*align(1) const V, @ptrCast(p)).*;
}

/// Stores `v` to `p`, unaligned; see `load`.
pub inline fn store(comptime V: type, p: [*]@typeInfo(V).vector.child, v: V) void {
    @as(*align(1) V, @ptrCast(p)).* = v;
}

/// The lanes of `v` summed as the vector kernels' reductions do: halves added
/// down to 128 bits (`vextractf128`, `vaddps`), then neighbors (`haddps`,
/// `haddpd`, NEON's `faddp`, SVE's `faddv`), so 4 lanes give `(v0+v1)+(v2+v3)`.
pub fn hsum(v: anytype) @typeInfo(@TypeOf(v)).vector.child {
    const F = @typeInfo(@TypeOf(v)).vector.child;
    const n = @typeInfo(@TypeOf(v)).vector.len;
    if (n == 1) return v[0];
    const h = n / 2;
    if (n * @sizeOf(F) > 16) {
        const lo = @shuffle(F, v, undefined, std.simd.iota(i32, h));
        const hi = @shuffle(F, v, undefined, std.simd.iota(i32, h) + @as(@Vector(h, i32), @splat(h)));
        return hsum(lo + hi);
    }
    const two: @Vector(h, i32) = @splat(2);
    const even = @shuffle(F, v, undefined, std.simd.iota(i32, h) * two);
    const odd = @shuffle(F, v, undefined, std.simd.iota(i32, h) * two + @as(@Vector(h, i32), @splat(1)));
    return hsum(even + odd);
}

/// The lanes of `v` summed by halves to the last lane, as AVX-512's
/// `_mm512_reduce_add_ps` and `_pd` do (gcc's and clang's), so 4 lanes give
/// `(v0+v2)+(v1+v3)`.
pub fn hsumFold(v: anytype) @typeInfo(@TypeOf(v)).vector.child {
    const F = @typeInfo(@TypeOf(v)).vector.child;
    const n = @typeInfo(@TypeOf(v)).vector.len;
    if (n == 1) return v[0];
    const h = n / 2;
    const lo = @shuffle(F, v, undefined, std.simd.iota(i32, h));
    const hi = @shuffle(F, v, undefined, std.simd.iota(i32, h) + @as(@Vector(h, i32), @splat(h)));
    return hsumFold(lo + hi);
}

/// Whether `@Vector` arithmetic in `F` keeps denormals here; false where kernels take
/// their scalar path (`comparisons.md`, 10): `f32` on ARMv7 NEON (Arm DDI 0601
/// 2026-09, AArch32 `FPSCR.FZ`) and AltiVec without VSX (`VSCR[NJ]`, as Linux v6.17
/// `arch/powerpc/kernel/process.c:1995` and QEMU v10.1.0 `target/ppc/cpu_init.c:4157` set it).
pub fn vectorIeee(comptime F: type) bool {
    return vectorIeeeOn(F, sysinfo.target.cpu);
}

/// `vectorIeee` for `cpu`, a core's of this target.
pub fn vectorIeeeOn(comptime F: type, comptime cpu: std.Target.Cpu) bool {
    if (F != f32 or sysinfo.target.abi.float() == .soft) return true;
    return switch (cpu.arch) {
        .arm, .armeb, .thumb, .thumbeb => !std.Target.arm.featureSetHas(cpu.features, .neon),
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => !std.Target.powerpc.featureSetHas(cpu.features, .altivec) or
            std.Target.powerpc.featureSetHas(cpu.features, .vsx),
        else => true,
    };
}

/// `x` (a float or an array of them) through an empty asm where AltiVec's `f32` lanes
/// flush, so LLVM cannot vectorize a loop through it; `x` itself elsewhere.
pub inline fn scalar(comptime T: type, x: T) T {
    if (comptime vectorIeee(f32) or !sysinfo.target.cpu.arch.isPowerpc()) return x;
    if (@typeInfo(T) == .array) {
        var y = x;
        inline for (&y) |*e| e.* = scalar(@typeInfo(T).array.child, e.*);
        return y;
    }
    if (T != f32) return x;
    return asm (""
        : [y] "=f" (-> f32),
        : [x] "0" (x),
    );
}

/// `x` rounded to `F` through memory where `F` computes on the x87, which
/// keeps a result wider until it is stored; `x` itself elsewhere. For a
/// value compared with one that was stored.
pub inline fn rounded(comptime F: type, x: F) F {
    if (comptime !onX87(F)) return x;
    var t: F = undefined;
    const v: *volatile F = &t;
    v.* = x;
    return v.*;
}

/// Whether `F` computes on the x87 here: x86 without SSE for `f32`, or
/// without SSE2 for `f64`.
pub fn onX87(comptime F: type) bool {
    return switch (F) {
        f32 => x87_f32,
        f64 => x87_f64,
        else => false,
    };
}

const x87_f32 = x87Without(.sse);
const x87_f64 = x87Without(.sse2);

fn x87Without(comptime unit: std.Target.x86.Feature) bool {
    const cpu = sysinfo.target.cpu;
    const has = std.Target.x86.featureSetHas;
    return cpu.arch.isX86() and has(cpu.features, .x87) and !has(cpu.features, unit);
}

/// Whether `n` elements of `size` bytes at stride `inc` fit in the address space
/// (`((n - 1) * |inc| + 1) * size` bytes, no overflow); the interfaces refuse the rest.
pub fn fits(n: isize, inc: isize, size: usize) bool {
    if (inc == std.math.minInt(isize)) return false;
    if (n <= 1) return true;
    return fitsMatrix(1, n, @intCast(@abs(inc)), size);
}

/// Whether a column-major `rows` by `cols` matrix with leading dimension `ld`
/// fits, as `fits` does for a vector: `(cols - 1) * ld + rows` elements.
pub fn fitsMatrix(rows: isize, cols: isize, ld: isize, size: usize) bool {
    if (rows <= 0 or cols <= 0) return true;
    const span = std.math.mul(isize, cols - 1, ld) catch return false;
    const elems = std.math.add(isize, span, rows) catch return false;
    _ = std.math.mul(isize, elems, std.math.cast(isize, size) orelse return false) catch return false;
    return true;
}

/// 64-bit `nrm2` for the portable and `@Vector` kernels: Blue's, as Reference-LAPACK 3.10's
/// `BLAS/SRC/dnrm2.f90` (Anderson, "Algorithm 978", ACM TOMS 44(1), 2017,
/// <https://doi.org/10.1145/3061665>); one scaled pass, no division; Inf and NaN propagate.
pub fn blue(comptime F: type, comptime w: comptime_int, n: isize, x: [*]F, inc_x: isize) F {
    // Fortran's `digits`, `minexponent` and `maxexponent`
    const digits = std.math.floatMantissaBits(F) + 1;
    const minexp = std.math.floatExponentMin(F) + 1;
    const maxexp = std.math.floatExponentMax(F) + 1;
    const tsml = std.math.ldexp(@as(F, 1), @divFloor(minexp, 2));
    const tbig = std.math.ldexp(@as(F, 1), @divFloor(maxexp - digits + 1, 2));
    const ssml = std.math.ldexp(@as(F, 1), -@divFloor(minexp - digits, 2));
    const sbig = std.math.ldexp(@as(F, 1), -@divFloor(maxexp + digits, 2));
    var notbig = true;
    var asml: F = 0;
    var amed: F = 0;
    var abig: F = 0;
    var i: isize = 0;
    while (i < n) : (i += 1) {
        inline for (0..w) |c| {
            const ax = @abs(at(F, x, i * inc_x * w + c).*);
            if (ax > tbig) {
                abig += (ax * sbig) * (ax * sbig);
                notbig = false;
            } else if (ax < tsml) {
                if (notbig) asml += (ax * ssml) * (ax * ssml);
            } else {
                amed += ax * ax;
            }
        }
    }
    // a NaN lands in `amed`
    const mixed = amed > 0 or amed != amed;
    if (abig > 0) {
        if (mixed) abig += (amed * sbig) * sbig;
        return @sqrt(abig) / sbig;
    }
    if (asml > 0) {
        if (!mixed) return @sqrt(asml) / ssml;
        const m = @sqrt(amed);
        const s = @sqrt(asml) / ssml;
        const lo = if (s > m) m else s;
        const hi = if (s > m) s else m;
        return @sqrt(hi * hi * (1 + (lo / hi) * (lo / hi)));
    }
    return @sqrt(amed);
}

/// `b / a` in parts by Smith's algorithm, as gfortran divides; no reciprocal, so a tiny
/// `a` still gives a finite quotient. Each step `rounded` on the x87.
pub inline fn cdiv(comptime F: type, br: F, bi: F, ar: F, ai: F) [2]F {
    const q = struct {
        inline fn f(x: F) F {
            return rounded(F, x);
        }
    }.f;
    if (@abs(ar) >= @abs(ai)) {
        const r = q(ai / ar);
        const d = q(ar + q(ai * r));
        return .{ q(q(br + q(bi * r)) / d), q(q(bi - q(br * r)) / d) };
    }
    const r = q(ar / ai);
    const d = q(ai + q(ar * r));
    return .{ q(q(q(br * r) + bi) / d), q(q(q(bi * r) - br) / d) };
}

/// `a * b` in parts, `std.math.Complex.mul`'s order, each step `rounded` on
/// the x87 as in `cdiv`; for products two routines promise the same bits of.
pub inline fn cmul(comptime F: type, ar: F, ai: F, br: F, bi: F) [2]F {
    const q = struct {
        inline fn f(x: F) F {
            return rounded(F, x);
        }
    }.f;
    return .{ q(q(ar * br) - q(ai * bi)), q(q(ai * br) + q(ar * bi)) };
}

/// What `extendedSum` adds up per float.
pub const Term = enum { abs, square };

/// x86_64 `asum.S`'s and `nrm2.S`'s x87 sums (`nrm2` before its root), for floats wider than
/// `f64`: 8 floats a trip into accumulator `float % 4`, the rest's float `c` into `c`,
/// folded `a3 + ((a2 + a0) + a1)` unfused; at `f80` the assembly's bits.
pub fn extendedSum(comptime F: type, comptime w: comptime_int, comptime term: Term, n: isize, x: [*]F, inc_x: isize) F {
    const per = 8 / w;
    var a: [4]F = @splat(0);
    var i: isize = 0;
    const main = n & -per;
    while (i < main) : (i += per) {
        inline for (0..8) |f| a[f % 4] += of(F, term, at(F, x, (i + f / w) * inc_x * w + f % w).*);
    }
    while (i < n) : (i += 1) {
        inline for (0..w) |c| a[c] += of(F, term, at(F, x, i * inc_x * w + c).*);
    }
    return a[3] + ((a[2] + a[0]) + a[1]);
}

inline fn of(comptime F: type, comptime term: Term, v: F) F {
    return switch (term) {
        .abs => @abs(v),
        .square => v * v,
    };
}
