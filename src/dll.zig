//! Run-time library loading: `Dll` per OS, file names, Darwin frameworks, and
//! `bind`, which fills a table of function pointers by symbol name.
//!
//! Lifetime: a handle whose pointers are kept for the rest of the process
//! (a static table, as `linalg.blas.apple_accel` keeps) is never closed;
//! closing it could unload code those pointers still reach.

const std = @import("std");
const sysinfo = @import("builtin");
const os_tag = sysinfo.target.os.tag;
const is_windows = os_tag == .windows;
const is_darwin = os_tag.isDarwin();
const DynLib = std.DynLib;
const win32 = @import("win32");
const Allocator = @import("zigstd.zig").Allocator;

/// Whether this target can load libraries at run time; elsewhere (freestanding,
/// UEFI, wasm) `Dll` opens nothing and finds nothing.
pub const supported = switch (os_tag) {
    .windows, .linux, .freebsd, .netbsd, .openbsd, .dragonfly, .illumos => true,
    else => is_darwin,
};

/// `dlopen` and `dlsym` are there: libc, not Linux's libc-free ELF loader.
const has_dl = supported and !is_windows and @hasField(@FieldType(DynLib, "inner"), "handle");

pub const WinDll = struct {
    const Self = @This();
    handle: win32.foundation.HINSTANCE,

    pub const Error = error{FileNotFound};

    /// Trusts the file. Malicious file will be able to execute arbitrary code.
    pub fn open(path: []const u8) Error!Self {
        var buf: [std.Io.Dir.max_path_bytes:0]u8 = undefined;
        const path_z = std.mem.printSentinel(&buf, "{s}", .{path}, 0) catch return error.FileNotFound;
        return openZ(path_z);
    }

    /// Trusts the file. Malicious file will be able to execute arbitrary code.
    pub fn openZ(path_c: [*:0]const u8) Error!Self {
        const handle = win32.kernel32.LoadLibraryA(path_c) orelse
            return error.FileNotFound;
        return .{ .handle = handle };
    }

    /// Trusts the file.
    pub fn close(self: *Self) void {
        _ = win32.kernel32.FreeLibrary(self.handle);
        self.* = undefined;
    }

    pub fn lookup(self: *Self, comptime T: type, name: [:0]const u8) ?T {
        const proc = win32.kernel32.GetProcAddress(self.handle, name.ptr) orelse
            return null;
        return @as(T, @ptrCast(proc));
    }
};

/// `Dll` where nothing can be loaded: every open fails.
pub const NoDll = struct {
    const Self = @This();

    pub const Error = error{FileNotFound};

    pub fn open(_: []const u8) Error!Self {
        return error.FileNotFound;
    }

    pub fn openZ(_: [*:0]const u8) Error!Self {
        return error.FileNotFound;
    }

    pub fn close(self: *Self) void {
        self.* = undefined;
    }

    pub fn lookup(_: *Self, comptime T: type, _: [:0]const u8) ?T {
        return null;
    }
};

pub const Dll = if (is_windows) WinDll else if (supported) DynLib else NoDll;

/// Returns a handle to `name` only if it is already loaded into this
/// process.
///
/// Do not call `.close()` on the result: freeing it would drop a
/// refcount this code never held, risking an early unload other code is still
/// depending on. Treat the returned handle as a permanent, borrowed reference
/// for the rest of the process's lifetime.
pub fn findLoaded(name: [:0]const u8) ?Dll {
    if (is_windows) {
        const handle = win32.kernel32.GetModuleHandleA(name.ptr) orelse
            return null;
        return .{ .handle = handle };
    }
    if (!has_dl) return null;
    const handle = std.c.dlopen(name.ptr, .{ .LAZY = true, .NOLOAD = true }) orelse
        return null;
    return .{ .inner = .{ .handle = handle } };
}

/// A library's file name on this OS, as a format of its base name.
const file_name_fmt = if (is_windows) "{s}.dll" else if (is_darwin) "lib{s}.dylib" else "lib{s}.so";

/// The file name of library `base` on this OS: `base.dll`, `libbase.dylib`
/// or `libbase.so`.
pub fn fileName(comptime base: []const u8) [:0]const u8 {
    return std.fmt.comptimePrint(file_name_fmt, .{base});
}

pub fn loadDll(alloc: Allocator, lib_name: []const u8) !Dll {
    const lib = try alloc.print(file_name_fmt, .{lib_name});
    defer alloc.free(lib);
    return try Dll.open(lib);
}

pub const DarwinFrameworkError = error{ FileNotFound, NameTooLong };

/// Opens system framework `name` (`/System/Library/Frameworks/<name>.framework/<name>`);
/// `FileNotFound` off Darwin. A framework is Apple's bundle of a shared library
/// with its headers and resources, on every Darwin OS. Since macOS 11 the file
/// is only in the dyld shared cache, which `dlopen` searches by that path.
pub fn openDarwinFramework(name: []const u8) DarwinFrameworkError!Dll {
    if (!has_dl or !is_darwin) return error.FileNotFound;
    var buf: [256:0]u8 = undefined;
    const path = std.mem.printSentinel(&buf, "/System/Library/Frameworks/{0s}.framework/{0s}", .{name}, 0) catch
        return error.NameTooLong;
    return Dll.openZ(path) catch error.FileNotFound;
}

/// What `bind` found: `fns`, `T` with the fields whose symbol resolved, and
/// which of each field's names did.
pub fn Bound(comptime T: type) type {
    return struct {
        const Self = @This();
        pub const Field = std.meta.FieldEnum(T);

        /// A field none of whose names resolved is null.
        fns: T,
        /// Per field, the index of the name that resolved in its list.
        matched: std.enums.EnumArray(Field, ?u8),

        /// The symbol name `field` was bound from.
        pub fn name(self: *const Self, comptime names: Names, comptime field: Field) ?[:0]const u8 {
            const list = comptime names(@tagName(field));
            const i = self.matched.get(field) orelse return null;
            inline for (list, 0..) |n, j| if (i == j) return n;
            unreachable;
        }

        /// How many fields were bound.
        pub fn count(self: *const Self) usize {
            var n: usize = 0;
            for (self.matched.values) |m| n += @intFromBool(m != null);
            return n;
        }
    };
}

/// A field's symbol names, tried in order.
pub const Names = fn (comptime field: []const u8) []const [:0]const u8;

/// The field's own name only.
pub fn exactName(comptime field: []const u8) []const [:0]const u8 {
    return &.{field ++ ""};
}

/// Fills each optional function-pointer field of `T` with the first of
/// `names(field)` that `lib` has; fields none resolve stay null.
pub fn bind(comptime T: type, lib: *Dll, comptime names: Names) Bound(T) {
    @setEvalBranchQuota(100_000);
    const info = @typeInfo(T).@"struct";
    var out: Bound(T) = .{ .fns = undefined, .matched = .initFill(null) };
    inline for (info.field_names, info.field_types) |f, F| {
        const P = switch (@typeInfo(F)) {
            .optional => |o| o.child,
            else => @compileError("bind: field " ++ f ++ " is not optional"),
        };
        if (@typeInfo(P) != .pointer or @typeInfo(@typeInfo(P).pointer.child) != .@"fn")
            @compileError("bind: field " ++ f ++ " is not a function pointer");
        const list = comptime names(f);
        if (list.len > std.math.maxInt(u8)) @compileError("bind: too many names for " ++ f);
        @field(out.fns, f) = null;
        inline for (list, 0..) |n, i| {
            if (lib.lookup(P, n)) |p| {
                @field(out.fns, f) = p;
                out.matched.set(@field(Bound(T).Field, f), i);
                break;
            }
        }
    }
    return out;
}

// ---- tests -------------------------------------------------------------------

const zigstd = @import("zigstd.zig");

/// A library every process here has loaded, and two of its functions.
const test_lib = switch (os_tag) {
    .windows => "kernel32.dll",
    .linux => "libc.so.6",
    else => "/usr/lib/libSystem.B.dylib",
};
const Cc = if (is_windows) std.lang.CallingConvention.winapi else std.lang.CallingConvention.c;
const TestTable = struct {
    pid: ?*const fn () callconv(Cc) u32,
    first: ?*const fn () callconv(Cc) u32,
    missing: ?*const fn () callconv(Cc) void,
};
const pid_fn = if (is_windows) "GetCurrentProcessId" else "getpid";
/// `first` lists two symbols that both exist; the first must win.
const first_fns = if (is_windows) .{ "GetCurrentThreadId", "GetCurrentProcessId" } else .{ "getuid", "getpid" };

fn testNames(comptime field: []const u8) []const [:0]const u8 {
    if (comptime std.mem.eql(u8, field, "pid")) return &.{ "rstd_no_such_symbol", pid_fn };
    if (comptime std.mem.eql(u8, field, "first")) return &.{ first_fns[0], first_fns[1] };
    return &.{ "rstd_no_such_symbol", "rstd_nor_this_one" };
}

fn testLibOk() bool {
    return switch (os_tag) {
        .windows => true,
        .linux => has_dl and sysinfo.target.abi.isGnu(),
        else => has_dl and is_darwin,
    };
}

test "bind fills found symbols in name order and leaves missing ones null" {
    if (comptime !testLibOk()) return error.SkipZigTest;
    var lib = try Dll.open(test_lib);
    defer lib.close();
    const b = bind(TestTable, &lib, testNames);
    try zigstd.testExpectTrue(b.fns.pid != null);
    try zigstd.testExpectEqual(@as(?u8, 1), b.matched.get(.pid));
    try zigstd.testExpectEqualStrings(pid_fn, b.name(testNames, .pid).?);
    try zigstd.testExpectEqual(@as(?u8, 0), b.matched.get(.first));
    try zigstd.testExpectTrue(b.fns.first == lib.lookup(@typeInfo(@FieldType(TestTable, "first")).optional.child, first_fns[0]).?);
    try zigstd.testExpectTrue(b.fns.missing == null);
    try zigstd.testExpectEqual(@as(?u8, null), b.matched.get(.missing));
    try zigstd.testExpectTrue(b.name(testNames, .missing) == null);
    try zigstd.testExpectEqual(@as(usize, 2), b.count());
    const pid: u32 = if (is_windows) win32.kernel32.GetCurrentProcessId() else @intCast(std.c.getpid());
    try zigstd.testExpectEqual(pid, b.fns.pid.?());
}

test "exact names, file names, and frameworks only on Darwin" {
    const want = if (is_windows) "m.dll" else if (is_darwin) "libm.dylib" else "libm.so";
    try zigstd.testExpectEqualStrings(want, fileName("m"));
    if (comptime !testLibOk()) return error.SkipZigTest;
    var lib = try Dll.open(test_lib);
    defer lib.close();
    const Exact = struct { getpid: ?*const fn () callconv(.c) c_int };
    if (!is_windows) try zigstd.testExpectTrue(bind(Exact, &lib, exactName).fns.getpid != null);
    if (is_darwin) {
        // kept open for the process, like any framework handle
        _ = try openDarwinFramework("CoreFoundation");
    } else {
        try zigstd.testExpectError(error.FileNotFound, openDarwinFramework("CoreFoundation"));
    }
    try zigstd.testExpectError(error.FileNotFound, Dll.open("rstd-no-such-library"));
}
