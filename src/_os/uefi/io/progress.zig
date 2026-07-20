const std = @import("std");
const uefi = std.os.uefi;
const Io = std.Io;

pub const progressParentFile = Io.failingProgressParentFile;
