//! `pty-tee [--conpty <conpty.dll>] <log> <command line>`: runs a command in a pseudoconsole the
//! size of this console, its output shown here and copied to `<log>` as written. `zig build
//! --watch` so sees a real console (its full progress tree, resized with the window) while another
//! reader follows the log. The command line passes through verbatim. Windows only; elsewhere
//! `script` does this.
//!
//! Windows' own pseudoconsole redraws the output; Windows Terminal's `conpty.dll` (1.22 on, as
//! WezTerm and VS Code ship it, with its `OpenConsole.exe`) passes the command's bytes through.
//! Without `--conpty`, the DLL search finds one beside pty-tee or on `PATH`.

const std = @import("std");
const sysinfo = @import("builtin");
const rstd = @import("rstd");

comptime {
    if (sysinfo.target.os.tag != .windows) @compileError("pty-tee is Windows only; use script elsewhere");
}

const win32 = rstd.win32;
const k32 = win32.kernel32;
const console = win32.system.console;
const threading = win32.system.threading;
const fs = win32.storage.file_system;
const HANDLE = win32.foundation.HANDLE;
const COORD = console.COORD;
const L = rstd.bytes.L;

/// The pseudoconsole calls, bound from a `conpty.dll` by name.
const Conpty = struct {
    CreatePseudoConsole: ?*const @TypeOf(k32.CreatePseudoConsole),
    ResizePseudoConsole: ?*const @TypeOf(k32.ResizePseudoConsole),
    ClosePseudoConsole: ?*const @TypeOf(k32.ClosePseudoConsole),

    const inbox: Conpty = .{
        .CreatePseudoConsole = &k32.CreatePseudoConsole,
        .ResizePseudoConsole = &k32.ResizePseudoConsole,
        .ClosePseudoConsole = &k32.ClosePseudoConsole,
    };

    /// `path`'s calls, or null when it does not load or lacks one.
    fn load(path: []const u8) ?Conpty {
        var lib = rstd.dll.Dll.open(path) catch return null;
        const bound = rstd.dll.bind(Conpty, &lib, rstd.dll.exactName);
        return if (bound.count() == @typeInfo(Conpty).@"struct".field_names.len) bound.fns else null;
    }
};

/// The size used when stdout is not a console.
const default_size: COORD = .{ .X = 120, .Y = 30 };

var out: ?HANDLE = null;
var log: HANDLE = undefined;
var pc: ?console.HPCON = null;
var conpty: Conpty = undefined;

const usage = "usage: pty-tee [--conpty <conpty.dll>] <log> <command line>\n";

pub fn main() void {
    const cmdline = k32.GetCommandLineW() orelse fail("no command line");
    var i = skipToken(cmdline, 0, null);
    var arg: [32768:0]u16 = undefined;
    var next = skipToken(cmdline, i, &arg);
    if (std.mem.eql(u16, std.mem.sliceTo(&arg, 0), L("--conpty"))) {
        i = next;
        next = skipToken(cmdline, i, &arg);
        if (next == i) usageExit();
        var path: [arg.len * 3]u8 = undefined;
        const n = std.unicode.utf16LeToUtf8(&path, std.mem.sliceTo(&arg, 0)) catch fail("the conpty path is not UTF-16");
        conpty = Conpty.load(path[0..n]) orelse fail("cannot load the conpty dll");
        i = next;
        next = skipToken(cmdline, i, &arg);
    } else conpty = Conpty.load("conpty.dll") orelse inbox: {
        std.debug.print("pty-tee: no conpty.dll beside pty-tee or on PATH; Windows' own pseudoconsole redraws the output, so the log is not the command's bytes\n", .{});
        break :inbox .inbox;
    };
    if (next == i or cmdline[next] == 0) usageExit();
    const log_path = arg;
    i = next;

    log = k32.CreateFileW(&log_path, fs.FILE_GENERIC_WRITE, .{ .READ = 1, .WRITE = 1, .DELETE = 1 }, null, fs.CREATE_ALWAYS, fs.FILE_ATTRIBUTE_NORMAL, null);
    if (log == win32.foundation.INVALID_HANDLE_VALUE) fail("cannot create the log");

    const in = k32.GetStdHandle(console.STD_INPUT_HANDLE);
    out = k32.GetStdHandle(console.STD_OUTPUT_HANDLE);
    var in_mode: console.CONSOLE_MODE = .{};
    var out_mode: console.CONSOLE_MODE = .{};
    const in_console = k32.GetConsoleMode(in, &in_mode) != 0;
    const out_console = k32.GetConsoleMode(out, &out_mode) != 0;
    const in_cp = k32.GetConsoleCP();
    const out_cp = k32.GetConsoleOutputCP();
    if (in_console) {
        // raw VT input: keys (Ctrl+C as 0x03) go to the pseudoconsole, which signals the command
        const cooked = bits(console.ENABLE_PROCESSED_INPUT) | bits(console.ENABLE_LINE_INPUT) | bits(console.ENABLE_ECHO_INPUT);
        _ = k32.SetConsoleMode(in, @bitCast((bits(in_mode) & ~cooked) | bits(console.ENABLE_VIRTUAL_TERMINAL_INPUT)));
        _ = k32.SetConsoleCP(win32.globalization.CP_UTF8);
    }
    if (out_console) {
        _ = k32.SetConsoleMode(out, @bitCast(bits(out_mode) | bits(console.ENABLE_PROCESSED_OUTPUT) |
            bits(console.ENABLE_VIRTUAL_TERMINAL_PROCESSING) | bits(console.DISABLE_NEWLINE_AUTO_RETURN)));
        _ = k32.SetConsoleOutputCP(win32.globalization.CP_UTF8);
    }
    // Ctrl+C ignored by `null` is inherited, so the command would ignore it too: clear any this
    // process inherited, and outlive Ctrl+C by a handler instead
    _ = k32.SetConsoleCtrlHandler(null, 0);
    _ = k32.SetConsoleCtrlHandler(&swallowCtrl, 1);

    var pty_in_r: ?HANDLE = null;
    var pty_in_w: ?HANDLE = null;
    var pty_out_r: ?HANDLE = null;
    var pty_out_w: ?HANDLE = null;
    if (k32.CreatePipe(&pty_in_r, &pty_in_w, null, 0) == 0 or k32.CreatePipe(&pty_out_r, &pty_out_w, null, 0) == 0) fail("cannot create pipes");
    // the cursor's place is asked of this console, rather than the screen cleared
    const flags: u32 = if (in_console and out_console) console.PSEUDOCONSOLE_INHERIT_CURSOR else 0;
    if (conpty.CreatePseudoConsole.?(consoleSize() orelse default_size, pty_in_r, pty_out_w, flags, &pc).failed) fail("cannot create a pseudoconsole");
    // the pseudoconsole holds its own copies; closing these lets the output pipe end
    _ = k32.CloseHandle(pty_in_r);
    _ = k32.CloseHandle(pty_out_w);

    var list_size: usize = 0;
    _ = k32.InitializeProcThreadAttributeList(null, 1, 0, &list_size);
    const list = std.heap.page_allocator.alloc(u8, list_size) catch fail("out of memory");
    if (k32.InitializeProcThreadAttributeList(list.ptr, 1, 0, &list_size) == 0 or
        k32.UpdateProcThreadAttribute(list.ptr, 0, threading.PROC_THREAD_ATTRIBUTE_PSEUDOCONSOLE, pc, @sizeOf(console.HPCON), null, null) == 0)
        fail("cannot attach the pseudoconsole");
    // no std handles of ours: a redirected one would be inherited in place of the pseudoconsole's
    var si: threading.STARTUPINFOEXW = .{ .StartupInfo = std.mem.zeroes(threading.STARTUPINFOW), .lpAttributeList = list.ptr };
    si.StartupInfo.cb = @sizeOf(threading.STARTUPINFOEXW);
    si.StartupInfo.dwFlags = threading.STARTF_USESTDHANDLES;
    var pi: threading.PROCESS_INFORMATION = undefined;
    if (k32.CreateProcessW(null, cmdline + i, null, null, 0, threading.EXTENDED_STARTUPINFO_PRESENT, null, null, &si.StartupInfo, &pi) == 0)
        fail("cannot start the command");
    _ = k32.CloseHandle(pi.hThread);

    const reader = std.Thread.spawn(.{}, pump, .{pty_out_r}) catch fail("cannot start a thread");
    {
        const t = std.Thread.spawn(.{}, forward, .{ in, pty_in_w }) catch fail("cannot start a thread");
        t.detach();
    }
    if (out_console) {
        const t = std.Thread.spawn(.{}, resize, .{}) catch fail("cannot start a thread");
        t.detach();
    }

    _ = k32.WaitForSingleObject(pi.hProcess, win32.system.windows_programming.INFINITE);
    var code: u32 = 1;
    _ = k32.GetExitCodeProcess(pi.hProcess, &code);
    // flushes the last output and ends the pipe
    conpty.ClosePseudoConsole.?(pc);
    reader.join();
    _ = k32.CloseHandle(log);
    if (in_console) {
        _ = k32.SetConsoleMode(in, in_mode);
        _ = k32.SetConsoleCP(in_cp);
    }
    if (out_console) {
        _ = k32.SetConsoleMode(out, out_mode);
        _ = k32.SetConsoleOutputCP(out_cp);
    }
    k32.ExitProcess(code);
}

/// A console mode as its bits, which zigwin32 names after input flags only.
fn bits(mode: console.CONSOLE_MODE) u32 {
    return @bitCast(mode);
}

/// The pseudoconsole's output, to stdout and the log, until it ends.
fn pump(from: ?HANDLE) void {
    var buf: [1 << 16]u8 = undefined;
    while (true) {
        var n: u32 = 0;
        if (k32.ReadFile(from, &buf, buf.len, &n, null) == 0 or n == 0) return;
        if (out) |o| writeAll(o, buf[0..n]);
        writeAll(log, buf[0..n]);
    }
}

/// Stdin, console or not, to the pseudoconsole until it ends.
fn forward(from: ?HANDLE, to: ?HANDLE) void {
    var buf: [4096]u8 = undefined;
    while (true) {
        var n: u32 = 0;
        if (k32.ReadFile(from, &buf, buf.len, &n, null) == 0) return;
        // a console read can return nothing; a pipe or file only at its end
        if (n == 0) {
            var mode: console.CONSOLE_MODE = .{};
            if (k32.GetConsoleMode(from, &mode) == 0) return;
            continue;
        }
        writeAll(to, buf[0..n]);
    }
}

/// Ctrl+C and Ctrl+Break reach the command as input; this process outlives them to clean up.
fn swallowCtrl(ctrl: u32) callconv(.winapi) i32 {
    return @intFromBool(ctrl == console.CTRL_C_EVENT or ctrl == console.CTRL_BREAK_EVENT);
}

/// The pseudoconsole kept this console's size.
fn resize() void {
    var last = consoleSize() orelse return;
    while (true) {
        k32.Sleep(200);
        const now = consoleSize() orelse continue;
        if (now.X != last.X or now.Y != last.Y) {
            _ = conpty.ResizePseudoConsole.?(pc, now);
            last = now;
        }
    }
}

fn consoleSize() ?COORD {
    var info: console.CONSOLE_SCREEN_BUFFER_INFO = undefined;
    if (k32.GetConsoleScreenBufferInfo(out, &info) == 0) return null;
    return .{ .X = info.srWindow.Right - info.srWindow.Left + 1, .Y = info.srWindow.Bottom - info.srWindow.Top + 1 };
}

fn writeAll(to: ?HANDLE, bytes: []const u8) void {
    var rest = bytes;
    while (rest.len > 0) {
        var n: u32 = 0;
        if (k32.WriteFile(to, rest.ptr, @intCast(rest.len), &n, null) == 0) return;
        rest = rest[n..];
    }
}

/// The index past the argument at `start` and the spaces after it, split as
/// `CommandLineToArgvW` does (quotes group; `\"` is a quote). With `into`, the
/// argument unquoted and 0-terminated.
fn skipToken(s: [*:0]const u16, start: usize, into: ?*[32768:0]u16) usize {
    var i = start;
    var len: usize = 0;
    var quoted = false;
    while (s[i] != 0 and (quoted or (s[i] != ' ' and s[i] != '\t'))) : (i += 1) {
        if (s[i] == '\\' and s[i + 1] == '"') {
            i += 1;
        } else if (s[i] == '"') {
            quoted = !quoted;
            continue;
        }
        if (into) |b| if (len < b.len) {
            b[len] = s[i];
            len += 1;
        };
    }
    if (into) |b| b[len] = 0;
    while (s[i] == ' ' or s[i] == '\t') i += 1;
    return i;
}

fn usageExit() noreturn {
    std.debug.print(usage, .{});
    k32.ExitProcess(2);
}

fn fail(what: []const u8) noreturn {
    std.debug.print("pty-tee: {s} (error {t})\n", .{ what, k32.GetLastError() });
    k32.ExitProcess(1);
}
