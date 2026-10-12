const std = @import("std");
const uefi = std.os.uefi;
const Io = std.Io;

pub const processExecutableOpen = Io.failingProcessExecutableOpen;
pub const processExecutablePath = Io.failingProcessExecutablePath;
pub const processCurrentPath = Io.failingProcessCurrentPath;
pub const processSetCurrentDir = Io.failingProcessSetCurrentDir;
pub const processSetCurrentPath = Io.failingProcessSetCurrentPath;
pub const processReplace = Io.failingProcessReplace;
pub const processSpawn = Io.failingProcessSpawn;
pub const childWait = Io.unreachableChildWait;
pub const childKill = Io.unreachableChildKill;
pub const inheritParentDir = Io.failingInheritParentDir;
pub const inheritParentFile = Io.failingInheritParentFile;
