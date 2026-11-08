//! MIPS UHI semihosting console and exit, for the `malta` rigs (`-semihosting`).

const UHI_exit = 1;
const UHI_write = 5;

inline fn uhi3(op: usize, a0: usize, a1: usize, a2: usize) void {
    asm volatile ("sdbbp 1"
        :
        : [op] "{$25}" (op),
          [a0] "{$4}" (a0),
          [a1] "{$5}" (a1),
          [a2] "{$6}" (a2),
        : .{ .memory = true, .r2 = true, .r3 = true });
}

pub fn writeByte(c: u8) void {
    var ch = c;
    uhi3(UHI_write, 1, @intFromPtr(&ch), 1);
}

pub fn exit(code: usize) noreturn {
    uhi3(UHI_exit, code, 0, 0);
    while (true) {}
}
