//! What the OS tells the non-x86 detectors about the CPU, as
//! `driver/others/dynamic_*.c` ask it: the auxiliary vector (`getauxval`,
//! FreeBSD's `elf_aux_info`), `/proc` and `/sys`, sysctl (BSDs, Darwin), and
//! Windows' `processorFeature` and registry. Elsewhere nothing is known,
//! and the callers take the generic kernels.

const std = @import("std");
const sysinfo = @import("builtin");

/// The auxiliary vector entries the detectors read.
pub const Tag = enum { hwcap, hwcap2, platform };

/// The OS's number for `tag`, if it has one: Linux and FreeBSD number them
/// differently (`HWCAP` is 16 and 25).
fn number(comptime tag: Tag) ?usize {
    const AT = std.elf.AT;
    return switch (tag) {
        .hwcap => if (@hasDecl(AT, "HWCAP")) AT.HWCAP else null,
        .hwcap2 => if (@hasDecl(AT, "HWCAP2")) AT.HWCAP2 else null,
        .platform => if (@hasDecl(AT, "PLATFORM")) AT.PLATFORM else null,
    };
}

/// Returns 0 when there is no auxiliary vector to read; either because
/// the OS has none, or because we are freestanding and nobody handed
/// us one. Callers treat 0 as "no capabilities known", which lands them on
/// the generic kernels.
pub fn get(comptime tag: Tag) usize {
    const n = comptime number(tag);
    if (comptime n == null) return 0;
    switch (comptime sysinfo.target.os.tag) {
        // with libc linked, std leaves `getauxval` to it
        .linux => return if (sysinfo.link_libc) std.c.getauxval(n.?) else std.os.linux.getauxval(n.?),
        // FreeBSD fills the same vector and answers through libc
        .freebsd => {
            if (comptime !sysinfo.link_libc) return 0;
            var value: c_ulong = 0;
            if (std.c.elf_aux_info(@intCast(n.?), &value, @sizeOf(c_ulong)) != 0) return 0;
            return value;
        },
        else => return 0,
    }
}

pub fn hwcap() usize {
    return get(.hwcap);
}

pub fn hwcap2() usize {
    return get(.hwcap2);
}

/// The NUL-terminated string `AT_PLATFORM` points at (Linux), if any.
pub fn platform() ?[]const u8 {
    const ptr = get(.platform);
    if (ptr == 0) return null;
    return std.mem.sliceTo(@as([*:0]const u8, @ptrFromInt(ptr)), 0);
}

/// Reads up to `buf.len` bytes of the file at `path` (`/proc`, `/sys`). Empty
/// when it cannot be opened, or off Linux. It goes straight to the raw
/// syscalls rather than through `std.fs`: no allocator, no `Io` instance, and
/// it runs at library init where neither is necessarily set up yet.
pub fn readFile(path: [*:0]const u8, buf: []u8) []const u8 {
    if (comptime sysinfo.target.os.tag != .linux) return buf[0..0];
    const linux = std.os.linux;
    const rc = linux.open(path, .{ .ACCMODE = .RDONLY }, 0);
    if (@as(isize, @bitCast(rc)) < 0) return buf[0..0];
    const fd: linux.fd_t = @intCast(rc);
    defer _ = linux.close(fd);
    var len: usize = 0;
    while (len < buf.len) {
        const got = linux.read(fd, buf[len..].ptr, buf.len - len);
        if (@as(isize, @bitCast(got)) <= 0) break;
        len += got;
    }
    return buf[0..len];
}

/// The value of the first `key: value` line of `/proc/cpuinfo` text whose key
/// starts with `key`, trimmed; null if there is none.
pub fn cpuinfoField(text: []const u8, key: []const u8) ?[]const u8 {
    var lines = std.mem.splitScalar(u8, text, '\n');
    while (lines.next()) |line| {
        if (!std.mem.startsWith(u8, line, key)) continue;
        const colon = std.mem.findScalar(u8, line, ':') orelse continue;
        return std.mem.trim(u8, line[colon + 1 ..], " \t\r");
    }
    return null;
}

/// A numeric sysctl (`hw.optional.arm.FEAT_SME` on Darwin), where the OS has
/// one and libc is linked. Null when it does not exist.
pub fn sysctlInt(comptime name: [:0]const u8) ?u64 {
    switch (comptime sysinfo.target.os.tag) {
        .driverkit, .ios, .maccatalyst, .macos, .tvos, .visionos, .watchos, .freebsd, .netbsd, .dragonfly => {},
        else => return null,
    }
    if (comptime !sysinfo.link_libc) return null;
    var value: u64 = 0;
    var len: usize = @sizeOf(u64);
    if (std.c.sysctlbyname(name.ptr, &value, &len, null, 0) != 0) return null;
    // narrower answers (int) fill the low bytes, as the reference's reads do
    return switch (len) {
        4 => @as(u32, @truncate(value)),
        else => value,
    };
}

/// Windows `IsProcessorFeaturePresent(feature)`, for the numbers
/// `std.os.windows.PF` stops short of (beyond 63 it is the system call).
pub fn processorFeature(feature: u32) bool {
    if (comptime sysinfo.target.os.tag != .windows) return false;
    const k32 = struct {
        extern "kernel32" fn IsProcessorFeaturePresent(f: u32) callconv(.winapi) i32;
    };
    return k32.IsProcessorFeaturePresent(feature) != 0;
}

/// A hex number as the sysfs files print it ("0x00000000411fd070"), the low 32
/// bits, as `strtoul` into an int gives.
pub fn parseHex(text: []const u8) ?u32 {
    var t = std.mem.trim(u8, text, " \t\r\n\x00");
    if (std.mem.startsWith(u8, t, "0x") or std.mem.startsWith(u8, t, "0X")) t = t[2..];
    var end: usize = 0;
    while (end < t.len and std.ascii.isHex(t[end])) end += 1;
    if (end == 0) return null;
    return @truncate(std.fmt.parseInt(u64, t[0..end], 16) catch return null);
}

/// How many CPUs `/sys/devices/system/cpu/possible` ("0-7", or "0") names.
pub fn cpuCountFromPossible(text: []const u8) usize {
    const t = std.mem.trim(u8, text, " \t\r\n\x00");
    if (!std.mem.startsWith(u8, t, "0-")) return 1;
    var end: usize = 2;
    while (end < t.len and std.ascii.isDigit(t[end])) end += 1;
    return (std.fmt.parseInt(usize, t[2..end], 10) catch 0) + 1;
}

/// A `REG_QWORD` under `HKEY_LOCAL_MACHINE`; null if absent or another type.
pub fn registryQword(comptime subkey: []const u8, comptime value: []const u8) ?u64 {
    if (comptime sysinfo.target.os.tag != .windows) return null;
    const advapi = struct {
        extern "advapi32" fn RegGetValueW(hkey: usize, subkey: [*:0]const u16, value: [*:0]const u16, flags: u32, kind: ?*u32, data: ?*anyopaque, size: ?*u32) callconv(.winapi) i32;
    };
    const HKEY_LOCAL_MACHINE: usize = 0x80000002;
    const RRF_RT_REG_QWORD: u32 = 0x40;
    // the 64-bit view from a 32-bit process too, not `WOW6432Node` (issue #264)
    const RRF_SUBKEY_WOW6464KEY: u32 = 0x10000;
    var raw: u64 = 0;
    var size: u32 = @sizeOf(u64);
    const rc = advapi.RegGetValueW(
        HKEY_LOCAL_MACHINE,
        std.unicode.utf8ToUtf16LeStringLiteral(subkey),
        std.unicode.utf8ToUtf16LeStringLiteral(value),
        RRF_RT_REG_QWORD | RRF_SUBKEY_WOW6464KEY,
        null,
        &raw,
        &size,
    );
    return if (rc == 0) raw else null;
}

test "auxv access is safe to call on any target" {
    // On non-Linux this is a compile-time-constant 0; on Linux without an
    // auxiliary vector it is a runtime 0. Either way it must not fault.
    _ = hwcap();
    _ = hwcap2();
    _ = platform();
}

test "readFile of a missing file is empty" {
    var buf: [16]u8 = undefined;
    try std.testing.expectEqual(@as(usize, 0), readFile("/definitely/not/here", &buf).len);
}

test "cpuinfoField takes the first matching line" {
    const text = "processor\t: 0\nCPU part\t: 0xc09\nCPU part\t: 0xc0f\nmodel name\t: Loongson-3A4000 @ 1.8GHz\n";
    try std.testing.expectEqualStrings("0xc09", cpuinfoField(text, "CPU part").?);
    try std.testing.expectEqualStrings("Loongson-3A4000 @ 1.8GHz", cpuinfoField(text, "model name").?);
    try std.testing.expect(cpuinfoField(text, "Hardware") == null);
}

test "sysctlInt and processorFeature answer on any target" {
    _ = sysctlInt("hw.ncpu");
    _ = processorFeature(0);
}

test "sysfs parsing" {
    try std.testing.expectEqual(@as(?u32, 0x411fd070), parseHex("0x00000000411fd070\n"));
    try std.testing.expectEqual(@as(?u32, 0x1000), parseHex("1000"));
    try std.testing.expectEqual(@as(?u32, null), parseHex("\n"));
    try std.testing.expectEqual(@as(usize, 8), cpuCountFromPossible("0-7\n"));
    try std.testing.expectEqual(@as(usize, 1), cpuCountFromPossible("0\n"));
    try std.testing.expectEqual(@as(usize, 1), cpuCountFromPossible(""));
}
