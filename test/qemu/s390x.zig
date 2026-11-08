//! s390x boot for `s390-ccw-virtio` (64-bit mode); console the SCLP ASCII console, exit a disabled-wait PSW (qemu powers off; the
//! summary carries the result). `_start` enables the vector facility and AFP registers (CR0 bits 45, 46) first.

comptime {
    asm (
        \\.section .text.boot
        \\.globl _start
        \\_start:
        \\  larl %r1, qemutest_cr0
        \\  stctg %c0, %c0, 0(%r1)
        \\  lg %r2, 0(%r1)
        \\  lgfi %r3, 0x60000
        \\  ogr %r2, %r3
        \\  stg %r2, 0(%r1)
        \\  lctlg %c0, %c0, 0(%r1)
        \\  larl %r15, __stack_top
        \\  aghi %r15, -160
        \\  larl %r1, __bss_start
        \\  larl %r2, __bss_end
        \\1:
        \\  clgr %r1, %r2
        \\  jhe 2f
        \\  mvghi 0(%r1), 0
        \\  la %r1, 8(%r1)
        \\  j 1b
        \\2:
        \\  brasl %r14, qemutest_main
        \\3:
        \\  j 3b
        \\.balign 8
        \\qemutest_cr0:
        \\  .quad 0
    );
}

/// SCCB header, then one event buffer (or the event mask body).
const Header = extern struct {
    length: u16,
    function_code: u8,
    control_mask: [3]u8,
    response_code: u16,
};

const cmd_write_event_mask: u32 = 0x00780005;
const cmd_write_event_data: u32 = 0x00760005;
const event_ascii_console: u8 = 0x1a;
const mask_ascii_console: u32 = 1 << (32 - event_ascii_console);

var sccb: [4096]u8 align(4096) = undefined;
var line: [256]u8 = undefined;
var line_len: usize = 0;
var mask_set = false;

fn servc(cmd: u32) void {
    asm volatile (".insn rre,0xb2200000,%[cmd],%[sccb]"
        :
        : [cmd] "d" (@as(u64, cmd)),
          [sccb] "d" (@intFromPtr(&sccb)),
        : .{ .memory = true });
}

fn setMask() void {
    @memset(&sccb, 0);
    const hdr: *Header = @ptrCast(&sccb);
    hdr.length = 8 + 4 + 16;
    const body: [*]u8 = @ptrFromInt(@intFromPtr(&sccb) + 8);
    const words: [*]u32 = @ptrCast(@alignCast(body + 4));
    body[2] = 0;
    body[3] = 4; // mask length in words is `u16` at offset 2: 4
    words[0] = mask_ascii_console; // cp_receive_mask
    words[1] = mask_ascii_console; // cp_send_mask
    words[2] = mask_ascii_console; // receive_mask
    words[3] = mask_ascii_console; // send_mask
    servc(cmd_write_event_mask);
    mask_set = true;
}

fn flush() void {
    if (line_len == 0) return;
    if (!mask_set) setMask();
    @memset(&sccb, 0);
    const hdr: *Header = @ptrCast(&sccb);
    const ev_len: u16 = @intCast(6 + line_len);
    hdr.length = 8 + ev_len;
    const ev: [*]u8 = @ptrFromInt(@intFromPtr(&sccb) + 8);
    ev[0] = @intCast(ev_len >> 8);
    ev[1] = @intCast(ev_len & 0xff);
    ev[2] = event_ascii_console;
    ev[3] = 0; // flags
    ev[4] = 0;
    ev[5] = 0;
    @memcpy(ev[6..][0..line_len], line[0..line_len]);
    servc(cmd_write_event_data);
    line_len = 0;
}

pub fn writeByte(c: u8) void {
    line[line_len] = c;
    line_len += 1;
    if (c == '\n' or line_len == line.len) flush();
}

pub fn exit(code: usize) noreturn {
    _ = code;
    flush();
    const psw = [2]u64{ 0x0002000180000000, 0 };
    asm volatile ("lpswe 0(%[p])"
        :
        : [p] "a" (&psw),
        : .{ .memory = true });
    while (true) {}
}
