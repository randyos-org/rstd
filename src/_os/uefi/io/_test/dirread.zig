//! `dirRead` on the UEFI rigs' boot volume, whose root holds only `EFI`
//! (`build/qemu/__root__.zig` lays it out).

const std = @import("std");
const uefi_io = @import("../__root__.zig");

test "UEFI dirRead lists the boot volume's root: one directory, EFI" {
    const io = uefi_io.io_inst;
    const dir = std.Io.Dir.cwd();
    var buf: [std.Io.Dir.Reader.min_buffer_len]u8 align(@alignOf(usize)) = undefined;
    var r: std.Io.Dir.Reader = .init(dir, &buf);
    var count: usize = 0;
    while (try r.next(io)) |e| {
        count += 1;
        try std.testing.expectEqualStrings("EFI", e.name);
        try std.testing.expectEqual(std.Io.File.Kind.directory, e.kind);
    }
    try std.testing.expectEqual(@as(usize, 1), count);
    // a second pass after a reset sees the same
    r.state = .reset;
    try std.testing.expect((try r.next(io)) != null);
}
