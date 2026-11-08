//! The calling thread's floating-point control state: the x87 control word
//! and MXCSR's control bits; FPCR; FPSCR's (arm, powerpc); `frm`; FCSR's
//! (mips, loongarch); FPC's; FSR's. Threads do not inherit it, so a parallel
//! region captures it on the caller and each worker task adopts it
//! (`threading.md`). Exception flags stay with each thread. Empty where the
//! target has no FPU or none of these.

const std = @import("std");
const sysinfo = @import("builtin");
const util = @import("util.zig");

const FpState = @This();

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
/// The self-hosted x86_64 backend takes no `(%reg)` operand; it gets `m`
/// operands, which LLVM misreads as inputs (issues.md).
const fp_mem_ops = sysinfo.zig_backend == .stage2_x86_64;
/// The e500's SPE, which replaces the FPU; newer Zig has no such feature
/// (upstream fe259973b6), Zig 0.17 still does.
const fp_spe = if (@hasField(std.Target.powerpc.Feature, "spe"))
    std.Target.powerpc.featureSetHas(fp_cpu.features, @field(std.Target.powerpc.Feature, "spe"))
else
    false;

const fp_reg: ?FpReg = switch (fp_cpu.arch) {
    .x86, .x86_64 => if (std.Target.x86.featureSetHas(fp_cpu.features, .sse)) .{ .Int = u32, .ctl = 0xffc0, .rounding = 0xe000, .toward_zero = 0xe000 } else null,
    .aarch64, .aarch64_be => if (std.Target.aarch64.featureSetHas(fp_cpu.features, .fp_armv8)) .{ .Int = u64, .ctl = std.math.maxInt(u64), .rounding = 0x1c0_0000, .toward_zero = 0x1c0_0000 } else null,
    .arm, .armeb, .thumb, .thumbeb => if (std.Target.arm.featureSetHas(fp_cpu.features, .vfp2sp)) .{ .Int = u32, .ctl = 0x07ff_ff00, .rounding = 0x1c0_0000, .toward_zero = 0x1c0_0000 } else null,
    .powerpc, .powerpcle, .powerpc64, .powerpc64le => if (fp_hard and std.Target.powerpc.featureSetHas(fp_cpu.features, .fpu) and !fp_spe) .{ .Int = u64, .ctl = 0xff, .rounding = 0x3, .toward_zero = 0x1 } else null,
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
    if (comptime fp_mem_ops) asm volatile ("fnstcw %[w]"
        : [w] "=m" (w),
    ) else asm volatile ("fnstcw (%[p])"
        :
        : [p] "r" (&w),
        : .{ .memory = true });
    return w;
}

fn setX87Word(w: u16) void {
    if (comptime fp_mem_ops) asm volatile ("fldcw %[w]"
        :
        : [w] "m" (w),
        : .{ .memory = true }) else asm volatile ("fldcw (%[p])"
        :
        : [p] "r" (&w),
        : .{ .memory = true });
}

fn readFp() fp_reg.?.Int {
    switch (fp_cpu.arch) {
        .x86, .x86_64 => {
            var v: u32 = undefined;
            if (comptime fp_mem_ops) asm volatile ("stmxcsr %[v]"
                : [v] "=m" (v),
            ) else asm volatile ("stmxcsr (%[p])"
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
        .x86, .x86_64 => if (comptime fp_mem_ops) asm volatile ("ldmxcsr %[v]"
            :
            : [v] "m" (v),
            : .{ .memory = true }) else asm volatile ("ldmxcsr (%[p])"
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
    return util.rounded(F, va.* + vb.*);
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
