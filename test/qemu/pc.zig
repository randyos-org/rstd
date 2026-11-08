//! The PC's port I/O, COM1 and power-off, for the x86 and x86_64 rigs.

pub inline fn outb(port: u16, value: u8) void {
    asm volatile ("outb %[value], %[port]"
        :
        : [value] "{al}" (value),
          [port] "{dx}" (port),
    );
}

pub inline fn outw(port: u16, value: u16) void {
    asm volatile ("outw %[value], %[port]"
        :
        : [value] "{ax}" (value),
          [port] "{dx}" (port),
    );
}

pub inline fn inb(port: u16) u8 {
    return asm volatile ("inb %[port], %[value]"
        : [value] "={al}" (-> u8),
        : [port] "{dx}" (port),
    );
}

const com1: u16 = 0x3f8;

/// `c` to COM1.
pub fn writeByte(c: u8) void {
    // the transmit holding register is empty
    while (inb(com1 + 5) & 0x20 == 0) {}
    outb(com1, c);
}

/// S5, at the PIIX4 PM control register and at the newer machines' one.
pub fn powerOff() noreturn {
    outw(0xb004, 0x2000);
    outw(0x604, 0x2000);
    while (true) asm volatile ("hlt");
}
