const std = @import("std");
const sysinfo = @import("builtin");
const Vtable = @import("Vtable.zig");

/// A vtable slot's name; too many fields for the default branch quota.
pub const Slot = blk: {
    @setEvalBranchQuota(10_000);
    break :blk std.meta.FieldEnum(Vtable);
};

/// The function type of a kernel slot.
pub fn SlotFn(comptime slot: Slot) type {
    return @typeInfo(@FieldType(Vtable, @tagName(slot))).optional.child;
}

/// Equivalent of C's `base[offset]` where `offset` may be negative -- BLAS
/// routines allow negative `inc_x`/`inc_y` (iterate a vector in reverse),
/// which does not map onto Zig's `usize`-indexed pointer/slice arithmetic
/// the way it does onto a C pointer (`x[ix]` with negative `ix` is just
/// `*(x + ix)` in C).
pub inline fn at(comptime T: type, base: [*]T, offset: i64) *T {
    const addr: usize = @intFromPtr(base);
    const byte_offset: i64 = offset * @as(i64, @sizeOf(T));
    const new_addr: usize = if (byte_offset >= 0)
        addr + @as(usize, @intCast(byte_offset))
    else
        addr - @as(usize, @intCast(-byte_offset));
    return @ptrFromInt(new_addr);
}

/// The symbol of `core`'s kernel for `slot` (`blaskern_daxpy_k_haswell`). The
/// per-core objects (`build/blaskern.zig`) export it, with hidden visibility so it
/// never leaves a shared library rstd is built into, and `kernels/<arch>.zig`'s
/// `install` looks it up by the same name.
pub fn kernelSymbol(comptime slot: Slot, comptime core: anytype) []const u8 {
    return "blaskern_" ++ @tagName(slot) ++ "_" ++ @tagName(core);
}

/// Loads a `V` from `p`, which is not assumed aligned to `V`'s natural
/// alignment -- e.g. a kernel walking a caller's buffer by an arbitrary
/// stride, where nothing guarantees the vector-width alignment `@Vector`
/// otherwise wants. (Lives here, not in `simd/`, because the per-core kernel
/// objects that use it are their own module rooted at `_kernel.zig` and
/// cannot reach outside `linalg/blas/`; see `kernels/__root__.zig`.)
pub inline fn load(comptime V: type, p: [*]const @typeInfo(V).vector.child) V {
    return @as(*align(1) const V, @ptrCast(p)).*;
}

/// Stores `v` to `p`, unaligned; see `load`.
pub inline fn store(comptime V: type, p: [*]@typeInfo(V).vector.child, v: V) void {
    @as(*align(1) V, @ptrCast(p)).* = v;
}

/// Whether `@Vector` arithmetic in `F` keeps denormals on this target. Not
/// `f32` on ARMv7 NEON, which flushes to zero whatever FPSCR says, nor on
/// AltiVec without VSX, which flushes under `VSCR[NJ]`, as Linux and QEMU
/// start it. Neither has `f64` lanes; LLVM splits those into scalar IEEE
/// operations, as it does every vector under a soft-float ABI. Kernels take
/// their scalar path where this is false (`comparisons.md`, 10). LLVM's
/// auto-vectorizer keeps off NEON's `f32` but not AltiVec's: the AltiVec
/// cores' objects are built without it (`cpu/_corespec.zig`), and the loops
/// built for the target pass their `f32` results through `scalar`.
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

/// `x` (a float, or an array of them), through an empty asm in a scalar
/// register where AltiVec's `f32` lanes flush: LLVM's vectorizers cannot
/// widen a loop or a group of operations through it, so a loop that passes
/// what it stores or accumulates through it keeps out of those lanes. `x`
/// itself everywhere else; ARMv7's vectorizers already keep off NEON's `f32`.
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

/// Whether a vector of `n` elements of `size` bytes at stride `inc` fits in
/// the address space: `((n - 1) * |inc| + 1) * size` bytes, with no overflow.
/// Where it does not, no such vector can exist; the interfaces refuse it so
/// no index computation overflows.
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

/// `nrm2`'s algorithm at 64 bits, for the portable and `@Vector` kernels:
/// Blue's, as Reference-LAPACK 3.10's `BLAS/SRC/dnrm2.f90` has it
/// (Anderson, "Algorithm 978: Safe Scaling in the Level 1 BLAS", ACM TOMS
/// 44(1), 2017): one pass summing the squares of small, medium and big
/// magnitudes apart, each scaled by a power of two so none overflows or
/// underflows; no division. Infinite for an infinity, NaN for a NaN.
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

/// `b / a` in parts, by Smith's algorithm (the larger part of `a` divided
/// out), as gfortran divides; for the triangular solves of levels 2 and 3.
/// No reciprocal is formed, so a tiny `a` whose reciprocal overflows still
/// gives a finite quotient wherever one exists. On the x87 every step is
/// `rounded`, so the bits do not hang on where the compiler spills.
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

/// x86_64 `asum.S`'s and `nrm2.S`'s x87 sums and their complex forms', before
/// `nrm2`'s square root, which every float wider than `f64` takes: eight
/// floats a trip into accumulators by float `% 4`, then each remaining
/// element's float `c` (of `w`) into accumulator `c`, folded
/// `a3 + ((a2 + a0) + a1)`, nothing fused. At `f80` it gives the assembly's
/// bits; `kernel/x86/`'s are in the same order.
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

// ---- floating-point control state ------------------------------------------------

/// The calling thread's floating-point control state: the x87 control word
/// and MXCSR's control bits; FPCR; FPSCR's (arm, powerpc); `frm`; FCSR's
/// (mips, loongarch); FPC's; FSR's. Threads do not inherit it, so a parallel
/// region captures it on the caller and each worker task adopts it
/// (`threading.md`). Exception flags stay with each thread. Empty where the
/// target has no FPU or none of these.
pub const FpState = struct {
    x87: if (fp_x87) u16 else void,
    ctl: if (fp_reg) |r| r.Int else void,

    /// Whether there is no state to carry.
    pub const empty = !fp_x87 and fp_reg == null;
    /// Whether Zig's own float code follows it: not under a soft-float ABI,
    /// whose arithmetic is compiler_rt's; only asm kernels' then.
    pub const governs = !empty and fp_hard;

    pub fn capture() FpState {
        var s: FpState = undefined;
        if (comptime fp_x87) s.x87 = x87Word();
        if (comptime fp_reg) |r| s.ctl = readFp() & r.ctl;
        return s;
    }

    /// Sets this thread's control bits to `s`'s, leaving its flags.
    pub fn apply(s: FpState) void {
        if (comptime fp_x87) setX87Word(s.x87);
        if (comptime fp_reg) |r| writeFp((readFp() & ~@as(r.Int, r.ctl | r.clear)) | s.ctl);
    }

    /// Puts `s` in force on this thread; returns the state it replaced if
    /// that differed, to `apply` when the task ends.
    pub fn adopt(s: FpState) ?FpState {
        const own = capture();
        if (own.eql(s)) return null;
        s.apply();
        return own;
    }

    pub fn eql(a: FpState, b: FpState) bool {
        return std.meta.eql(a, b);
    }

    /// `s` rounding toward zero, and flushing denormals on x86's SSE and
    /// Arm's; for tests.
    pub fn towardZero(s: FpState) FpState {
        var t = s;
        if (comptime fp_x87) t.x87 |= 0xc00;
        if (comptime fp_reg) |r| t.ctl = (t.ctl & ~@as(r.Int, r.rounding)) | r.toward_zero;
        return t;
    }
};

/// A register of floating-point control and status bits.
const FpReg = struct {
    Int: type,
    /// Its control bits.
    ctl: comptime_int,
    /// Status bits `apply` clears: a cause bit set with its trap enabled traps.
    clear: comptime_int = 0,
    /// The rounding field (and flush to zero), and its value toward zero.
    rounding: comptime_int,
    toward_zero: comptime_int,
};

const fp_cpu = sysinfo.target.cpu;
const fp_hard = sysinfo.target.abi.float() == .hard;
const fp_x87 = fp_cpu.arch.isX86() and std.Target.x86.featureSetHas(fp_cpu.features, .x87);

const fp_reg: ?FpReg = switch (fp_cpu.arch) {
    .x86, .x86_64 => if (std.Target.x86.featureSetHas(fp_cpu.features, .sse)) .{ .Int = u32, .ctl = 0xffc0, .rounding = 0xe000, .toward_zero = 0xe000 } else null,
    .aarch64, .aarch64_be => if (std.Target.aarch64.featureSetHas(fp_cpu.features, .fp_armv8)) .{ .Int = u64, .ctl = std.math.maxInt(u64), .rounding = 0x1c0_0000, .toward_zero = 0x1c0_0000 } else null,
    .arm, .armeb, .thumb, .thumbeb => if (std.Target.arm.featureSetHas(fp_cpu.features, .vfp2sp)) .{ .Int = u32, .ctl = 0x07ff_ff00, .rounding = 0x1c0_0000, .toward_zero = 0x1c0_0000 } else null,
    .powerpc, .powerpcle, .powerpc64, .powerpc64le => if (fp_hard and std.Target.powerpc.featureSetHas(fp_cpu.features, .fpu) and !std.Target.powerpc.featureSetHas(fp_cpu.features, .spe)) .{ .Int = u64, .ctl = 0xff, .rounding = 0x3, .toward_zero = 0x1 } else null,
    .riscv32, .riscv32be, .riscv64, .riscv64be => if (std.Target.riscv.featureSetHas(fp_cpu.features, .f)) .{ .Int = usize, .ctl = 0x7, .rounding = 0x7, .toward_zero = 0x1 } else null,
    .mips, .mipsel, .mips64, .mips64el => if (fp_hard) .{ .Int = u32, .ctl = 0x0100_0f83, .clear = 0x0003_f000, .rounding = 0x3, .toward_zero = 0x1 } else null,
    .loongarch32, .loongarch64 => if (std.Target.loongarch.featureSetHas(fp_cpu.features, .f)) .{ .Int = u32, .ctl = 0x31f, .clear = 0x1f00_0000, .rounding = 0x300, .toward_zero = 0x100 } else null,
    .s390x => if (fp_hard) .{ .Int = u32, .ctl = 0xff00_0077, .rounding = 0x7, .toward_zero = 0x1 } else null,
    .sparc => if (fp_hard) .{ .Int = u32, .ctl = 0xcfc0_0000, .clear = 0x1f, .rounding = 0xc000_0000, .toward_zero = 0x4000_0000 } else null,
    .sparc64 => if (fp_hard) .{ .Int = u64, .ctl = 0xcfc0_0000, .clear = 0x1f, .rounding = 0xc000_0000, .toward_zero = 0x4000_0000 } else null,
    // no-op: none that rstd starts threads on
    else => null,
};

fn x87Word() u16 {
    var w: u16 = undefined;
    asm volatile ("fnstcw (%[p])"
        :
        : [p] "r" (&w),
        : .{ .memory = true });
    return w;
}

fn setX87Word(w: u16) void {
    asm volatile ("fldcw (%[p])"
        :
        : [p] "r" (&w),
        : .{ .memory = true });
}

fn readFp() fp_reg.?.Int {
    switch (fp_cpu.arch) {
        .x86, .x86_64 => {
            var v: u32 = undefined;
            asm volatile ("stmxcsr (%[p])"
                :
                : [p] "r" (&v),
                : .{ .memory = true });
            return v;
        },
        .aarch64, .aarch64_be => return asm volatile ("mrs %[r], fpcr"
            : [r] "=r" (-> u64),
        ),
        .arm, .armeb, .thumb, .thumbeb => return asm volatile ("vmrs %[r], fpscr"
            : [r] "=r" (-> u32),
        ),
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => return @bitCast(asm volatile ("mffs %[f]"
            : [f] "=f" (-> f64),
        )),
        .riscv32, .riscv32be, .riscv64, .riscv64be => return asm volatile ("frrm %[r]"
            : [r] "=r" (-> usize),
        ),
        .mips, .mipsel, .mips64, .mips64el => return asm volatile ("cfc1 %[r], $31"
            : [r] "=r" (-> u32),
        ),
        .loongarch32, .loongarch64 => return asm volatile ("movfcsr2gr %[r], $fcsr0"
            : [r] "=r" (-> u32),
        ),
        .s390x => return asm volatile ("efpc %[r]"
            : [r] "=r" (-> u32),
        ),
        .sparc, .sparc64 => {
            var v: fp_reg.?.Int = undefined;
            asm volatile ((if (fp_cpu.arch == .sparc64) "stx" else "st") ++ " %%fsr, [%[p]]"
                :
                : [p] "r" (&v),
                : .{ .memory = true });
            return v;
        },
        else => comptime unreachable,
    }
}

fn writeFp(v: fp_reg.?.Int) void {
    switch (fp_cpu.arch) {
        .x86, .x86_64 => asm volatile ("ldmxcsr (%[p])"
            :
            : [p] "r" (&v),
            : .{ .memory = true }),
        .aarch64, .aarch64_be => asm volatile ("msr fpcr, %[r]"
            :
            : [r] "r" (v),
            : .{ .memory = true }),
        .arm, .armeb, .thumb, .thumbeb => asm volatile ("vmsr fpscr, %[r]"
            :
            : [r] "r" (v),
            : .{ .memory = true }),
        // fields 6 and 7: the enables, NI and RN
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => asm volatile ("mtfsf 3, %[f]"
            :
            : [f] "f" (@as(f64, @bitCast(v))),
            : .{ .memory = true }),
        .riscv32, .riscv32be, .riscv64, .riscv64be => asm volatile ("fsrm %[r]"
            :
            : [r] "r" (v),
            : .{ .memory = true }),
        .mips, .mipsel, .mips64, .mips64el => asm volatile ("ctc1 %[r], $31"
            :
            : [r] "r" (v),
            : .{ .memory = true }),
        .loongarch32, .loongarch64 => asm volatile ("movgr2fcsr $fcsr0, %[r]"
            :
            : [r] "r" (v),
            : .{ .memory = true }),
        .s390x => asm volatile ("sfpc %[r]"
            :
            : [r] "r" (v),
            : .{ .memory = true }),
        .sparc, .sparc64 => asm volatile ((if (fp_cpu.arch == .sparc64) "ldx" else "ld") ++ " [%[p]], %%fsr"
            :
            : [p] "r" (&v),
            : .{ .memory = true }),
        else => comptime unreachable,
    }
}

/// `1 + 3/4 ulp`, rounded by this thread's mode: up to nearest, down toward zero.
fn onePlus(comptime F: type) F {
    var a: F = 1;
    var b: F = std.math.floatEps(F) * 0.75;
    const va: *volatile F = &a;
    const vb: *volatile F = &b;
    return rounded(F, va.* + vb.*);
}

/// Whether `f32` and `f64` round `onePlus` down.
fn roundsDown() [2]bool {
    return .{ onePlus(f32) == 1, onePlus(f64) == 1 };
}

test "FpState: apply sets what capture reads, the rounding follows, and the state comes back" {
    if (comptime FpState.empty) return error.SkipZigTest;
    const own: FpState = .capture();
    defer own.apply();
    const down = [2]bool{ FpState.governs, FpState.governs };
    try std.testing.expectEqual([2]bool{ false, false }, roundsDown());
    const rz = own.towardZero();
    try std.testing.expect(!rz.eql(own));
    try std.testing.expect(rz.adopt().?.eql(own));
    try std.testing.expect(FpState.capture().eql(rz));
    try std.testing.expect(rz.adopt() == null);
    try std.testing.expectEqual(down, roundsDown());
    own.apply();
    try std.testing.expect(FpState.capture().eql(own));
    try std.testing.expectEqual([2]bool{ false, false }, roundsDown());
}
