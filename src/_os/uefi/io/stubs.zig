//! The `Io` entries UEFI has no backing for: networking, processes and
//! progress reporting.

const std = @import("std");
const Io = std.Io;

pub const netListenIp = Io.failingNetListenIp;
pub const netBindIp = Io.failingNetBindIp;
pub const netConnectIp = Io.failingNetConnectIp;
pub const netListenUnix = Io.failingNetListenUnix;
pub const netConnectUnix = Io.failingNetConnectUnix;
pub const netSocketCreatePair = Io.failingNetSocketCreatePair;
pub const netWriteFile = Io.failingNetWriteFile;
pub const netClose = Io.unreachableNetClose;
pub const netShutdown = Io.failingNetShutdown;
pub const netInterfaceNameResolve = Io.failingNetInterfaceNameResolve;
pub const netInterfaceName = Io.unreachableNetInterfaceName;
pub const netLookup = Io.failingNetLookup;

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

pub const progressParentFile = Io.failingProgressParentFile;
