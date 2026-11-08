//! The Open Firmware client interface, for the rigs OpenBIOS or VOF starts:
//! ppc32 (32-bit cells) and sparc64 (64-bit).

/// The interface at `entry.*` (the address `_start` stored), in `Cell`s.
pub fn Client(comptime Cell: type, comptime Ret: type, comptime entry: *const usize) type {
    return struct {
        /// One call: the service, then its arguments and room for its results.
        pub fn call(cells: []Cell) void {
            const f: *const fn ([*]Cell) callconv(.c) Ret = @ptrFromInt(entry.*);
            _ = f(cells.ptr);
        }

        pub fn serviceName(comptime s: []const u8) Cell {
            return @intFromPtr((s ++ "\x00").ptr);
        }

        /// `c` to the device `handle`.
        pub fn write(handle: u32, c: u8) void {
            var ch = c;
            var cells = [_]Cell{ serviceName("write"), 3, 1, handle, @intFromPtr(&ch), 1, 0 };
            call(&cells);
        }

        /// Resets the machine, which `-no-reboot` turns into qemu exiting.
        pub fn resetAll() void {
            var reset = [_]Cell{ serviceName("interpret"), 1, 1, serviceName("reset-all"), 0 };
            call(&reset);
        }
    };
}
