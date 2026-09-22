//! Converts the inline-assembly microkernels of OpenBLAS's C `axpy` kernels
//! (`<p>axpy_microk_<core>-2.c`) into Zig functions.
//!
//!     OPENBLAS=<checkout> zig run tools/blas_asm/microk.zig
//!
//! A microkernel is a C function around one `__asm__` statement, or two
//! chosen by the length. The assembly is kept as written; only its spelling
//! changes, from GCC's to Zig's: numbered operands `%2` become named ones
//! `%[x]`, and the clobbers are declared as Zig spells them, at the width of
//! the vector registers the text actually uses. The `#if !defined(CONJ)`
//! inside the strings is resolved by the preprocessor, once for each form.
//!
//! The driver around a microkernel (`daxpy.c` and its siblings) is not
//! converted here; `chain.zig` reads which core uses which block, and the
//! architecture's `kernel.zig` composes it with the driver in
//! `axpy/kernel.zig`.

const std = @import("std");
const Io = std.Io;
const map = @import("kernel_map.zig");

const Precision = struct { letter: u8, float: []const u8, complex: bool };

fn precisionOf(letter: u8) Precision {
    return switch (letter) {
        's' => .{ .letter = 's', .float = "f32", .complex = false },
        'd' => .{ .letter = 'd', .float = "f64", .complex = false },
        'c' => .{ .letter = 'c', .float = "f32", .complex = true },
        else => .{ .letter = 'z', .float = "f64", .complex = true },
    };
}

fn put(b: *std.ArrayList(u8), comptime fmt: []const u8, args: anytype) !void {
    try b.print(map.arena, fmt, args);
}

/// Extra preprocessor arguments of the target being converted.
var cpp_extra: []const []const u8 = &.{};

/// Whether the target is POWER, whose blocks use VSX registers by number.
var is_power = false;

fn runCpp(source: []const u8, conj: bool) ![]const u8 {
    var argv: std.ArrayList([]const u8) = .empty;
    try argv.appendSlice(map.arena, &.{ "zig", "cc", "-E", "-P", "-x", "c", "-w" });
    try argv.appendSlice(map.arena, cpp_extra);
    if (conj) try argv.append(map.arena, "-DCONJ");
    try argv.appendSlice(map.arena, &.{ source, "-o", "-" });
    const r = std.process.run(map.arena, map.io, .{ .argv = argv.items }) catch |e|
        map.die("cannot run the preprocessor: {t}", .{e});
    return r.stdout;
}

/// The text of the string literals in `s`, joined, with C escapes resolved.
fn literals(s: []const u8) ![]const u8 {
    var out: std.ArrayList(u8) = .empty;
    var i: usize = 0;
    while (i < s.len) : (i += 1) {
        if (s[i] != '"') continue;
        i += 1;
        while (i < s.len and s[i] != '"') : (i += 1) {
            if (s[i] == '\\' and i + 1 < s.len) {
                i += 1;
                try out.append(map.arena, switch (s[i]) {
                    'n' => '\n',
                    't' => '\t',
                    else => s[i],
                });
            } else try out.append(map.arena, s[i]);
        }
    }
    return out.items;
}

/// LLVM writes a two-letter constraint such as GCC's `wa` (any VSX register)
/// with a `^` in front.
fn spellConstraint(c: []const u8) ![]const u8 {
    if (!is_power) return c;
    const at = std.mem.indexOfScalar(u8, c, 'w') orelse return c;
    if (at + 1 >= c.len) return c;
    return try std.fmt.allocPrint(map.arena, "{s}^{s}", .{ c[0..at], c[at..] });
}

/// One `"constraint" (expression)` operand.
const Operand = struct { constraint: []const u8, expr: []const u8, name: []const u8, dropped: bool = false };

fn parseOperands(s: []const u8, first_index: usize) ![]const Operand {
    var out: std.ArrayList(Operand) = .empty;
    var depth: usize = 0;
    var start: usize = 0;
    var i: usize = 0;
    while (i <= s.len) : (i += 1) {
        const c: u8 = if (i < s.len) s[i] else ',';
        if (c == '(') depth += 1;
        if (c == ')') depth -= 1;
        if (c == ',' and depth == 0) {
            const part = std.mem.trim(u8, s[start..i], " \t\r\n");
            start = i + 1;
            if (part.len == 0) continue;
            const q1 = std.mem.indexOfScalar(u8, part, '"').?;
            const q2 = std.mem.indexOfScalarPos(u8, part, q1 + 1, '"').?;
            const open = std.mem.indexOfScalarPos(u8, part, q2, '(').?;
            const close = std.mem.lastIndexOfScalar(u8, part, ')').?;
            const expr = std.mem.trim(u8, part[open + 1 .. close], " ");
            // a written name is the operand's; otherwise a plain variable
            // names it, and anything else is numbered
            const plain = !std.ascii.isDigit(expr[0]) and for (expr) |ch| {
                if (!std.ascii.isAlphanumeric(ch) and ch != '_') break false;
            } else true;
            const name = if (part[0] == '[')
                part[1..std.mem.indexOfScalar(u8, part, ']').?]
            else if (plain) expr else try std.fmt.allocPrint(map.arena, "op{d}", .{first_index + out.items.len});
            // a memory operand the C describes by an array cast has no Zig
            // spelling; the `memory` clobber says the same
            const drop = std.mem.indexOf(u8, expr, "*(") != null and std.mem.indexOfScalar(u8, part[q1 + 1 .. q2], 'm') != null;
            try out.append(map.arena, .{ .constraint = try spellConstraint(part[q1 + 1 .. q2]), .expr = expr, .name = name, .dropped = drop });
        }
    }
    return out.items;
}

/// The operand sections of an `asm` statement: what follows its template,
/// split at the top-level colons.
fn sections(s: []const u8) [3][]const u8 {
    var out: [3][]const u8 = .{ "", "", "" };
    var idx: usize = 0;
    var start: usize = 0;
    var i: usize = 0;
    var in_str = false;
    var depth: usize = 0;
    while (i < s.len) : (i += 1) {
        const c = s[i];
        if (c == '"' and (i == 0 or s[i - 1] != '\\')) in_str = !in_str;
        if (in_str) continue;
        if (c == '(') depth += 1;
        if (c == ')') depth -= 1;
        if (c == ':' and depth == 0) {
            if (idx > 0 and idx <= 3) out[idx - 1] = s[start..i];
            idx += 1;
            start = i + 1;
        }
    }
    if (idx > 0 and idx <= 3) out[idx - 1] = s[start..];
    return out;
}

/// The Zig type of a C parameter.
fn zigType(c_type: []const u8, float: []const u8, name: []const u8) ![]const u8 {
    if (std.mem.eql(u8, c_type, "BLASLONG") or std.mem.eql(u8, c_type, "long")) return "isize";
    if (std.mem.eql(u8, c_type, "double") or std.mem.eql(u8, c_type, "float")) return float;
    if (std.mem.eql(u8, c_type, "FLOAT *")) {
        const c = if (std.mem.eql(u8, name, "alpha")) "const " else "";
        return try std.fmt.allocPrint(map.arena, "[*]{s}{s}", .{ c, float });
    }
    map.die("unknown parameter type '{s}'", .{c_type});
}

/// Whether an operand `"+r" (name)` or `"+&r"(name)` writes the variable back.
fn isWritten(text: []const u8, name: []const u8) !bool {
    const needle = try std.fmt.allocPrint(map.arena, "({s})", .{name});
    var at: usize = 0;
    while (std.mem.indexOfPos(u8, text, at, needle)) |p| {
        const before = text[if (p > 12) p - 12 else 0..p];
        if (std.mem.indexOf(u8, before, "\"+") != null) return true;
        at = p + 1;
    }
    return false;
}

/// What the target calls the condition codes in a clobber list.
var cc_name: []const u8 = "cc";

/// A C expression in Zig: `*p` is `p[0]`.
fn zigExpr(e: []const u8) ![]const u8 {
    if (e.len > 1 and e[0] == '*') return try std.fmt.allocPrint(map.arena, "{s}[0]", .{e[1..]});
    return e;
}

/// The output an input constraint is tied to, when it is a number.
fn tieOf(constraint: []const u8) ?usize {
    return std.fmt.parseInt(usize, constraint, 10) catch null;
}

/// Whether some input is tied to output `k`.
fn isTied(ins: []const Operand, k: usize) bool {
    for (ins) |o| {
        if (tieOf(o.constraint) == k) return true;
    }
    return false;
}

fn allDropped(ops: []const Operand) bool {
    for (ops) |o| if (!o.dropped) return false;
    return true;
}

/// One `asm volatile` statement.
fn emitAsm(w: *std.ArrayList(u8), indent: []const u8, stmt: []const u8, locals: []const []const u8, name: []const u8) !void {
    const secs = sections(stmt);
    const colon = std.mem.indexOf(u8, stmt, secs[0]).? - 1;
    const template = try literals(stmt[0..colon]);
    // text outside the string literals is an unexpanded macro
    {
        var in_lit = false;
        for (stmt[0..colon], 0..) |c, k| {
            if (c == '"' and (k == 0 or stmt[k - 1] != '\\')) in_lit = !in_lit;
            if (!in_lit and c != '"' and !std.ascii.isWhitespace(c)) map.die("unexpanded macro text in the template of {s}: {s}", .{ name, stmt[k -| 40..k] });
        }
    }
    const outs = try parseOperands(secs[0], 0);
    const ins = try parseOperands(secs[1], outs.len);
    var all: std.ArrayList(Operand) = .empty;
    try all.appendSlice(map.arena, outs);
    try all.appendSlice(map.arena, ins);

    // `%N` to `%[name]`; `%%` is the register prefix and stays
    var text_out: std.ArrayList(u8) = .empty;
    var i: usize = 0;
    while (i < template.len) : (i += 1) {
        if (template[i] == '%' and i + 1 < template.len) {
            if (template[i + 1] == '%') {
                try text_out.appendSlice(map.arena, "%%");
                i += 1;
            } else if (template[i + 1] == '=') {
                // unique per expansion; Zig spells it the same
                try text_out.appendSlice(map.arena, "%=");
                i += 1;
            } else if (std.ascii.isDigit(template[i + 1]) or (template[i + 1] == 'x' and i + 2 < template.len and std.ascii.isDigit(template[i + 2]))) {
                // `%N`, and `%xN` (the VSX name of a register operand)
                const modifier = template[i + 1] == 'x';
                var j = i + 1 + @as(usize, @intFromBool(modifier));
                var n: usize = 0;
                while (j < template.len and std.ascii.isDigit(template[j])) : (j += 1) n = n * 10 + (template[j] - '0');
                // a dropped memory operand is only ever named in a comment
                if (all.items[n].dropped) {
                    try text_out.appendSlice(map.arena, "<mem>");
                } else try text_out.print(map.arena, "%[{s}{s}]", .{ all.items[n].name, if (modifier) ":x" else "" });
                i = j - 1;
            } else if (template[i + 1] == '[') {
                // a named operand, spelled as Zig does
                try text_out.append(map.arena, '%');
            } else map.die("unhandled operand form in {s}", .{name});
        } else try text_out.append(map.arena, template[i]);
    }

    // Zig takes 15 outputs in practice (16 crashes AstGen). The scratch vector
    // registers past that go in as inputs, which the asm overwrites and which
    // are dead after it; the outputs are then early-clobber so that none
    // shares a register. An output an input is tied to (`"16"`) stays.
    const moved_flag = try map.arena.alloc(bool, outs.len);
    @memset(moved_flag, false);
    var kept: usize = 0;
    for (outs) |o| {
        if (!o.dropped) kept += 1;
    }
    const spill = kept > 15;
    {
        var need: usize = if (spill) kept - 15 else 0;
        var k = outs.len;
        while (k > 0 and need > 0) {
            k -= 1;
            if (outs[k].dropped or isTied(ins, k)) continue;
            moved_flag[k] = true;
            need -= 1;
        }
    }
    // where each output lands among those emitted, which is what a tie names
    const position = try map.arena.alloc(usize, outs.len);
    {
        var p: usize = 0;
        for (outs, 0..) |o, k| {
            position[k] = p;
            if (!o.dropped and !moved_flag[k]) p += 1;
        }
    }
    // Each gets its own value, or the compiler shares one register among
    // the inputs that are all undefined.
    var distinct: usize = 1;
    for (outs, 0..) |o, k| {
        if (!moved_flag[k]) continue;
        try put(w, "{s}{s} = @splat({d});\n", .{ indent, o.expr, distinct });
        distinct += 1;
    }
    try put(w, "{s}asm volatile (\n", .{indent});
    var tl = std.mem.splitScalar(u8, text_out.items, '\n');
    while (tl.next()) |raw| {
        // tabs are not allowed in a Zig multiline string; the alignment is only cosmetic
        var line: std.ArrayList(u8) = .empty;
        var col: usize = 0;
        for (std.mem.trimEnd(u8, raw, " \t\r")) |c| {
            if (c == '\t') {
                const n = 8 - col % 8;
                try line.appendNTimes(map.arena, ' ', n);
                col += n;
            } else {
                try line.append(map.arena, c);
                col += 1;
            }
        }
        if (std.mem.trim(u8, line.items, " ").len == 0) continue;
        try put(w, "{s}    \\\\{s}\n", .{ indent, line.items });
    }

    var sep: []const u8 = ":";
    var moved: std.ArrayList(Operand) = .empty;
    for (outs, 0..) |o, k| {
        if (o.dropped) continue;
        if (moved_flag[k]) {
            try moved.append(map.arena, .{ .constraint = o.constraint[1..], .expr = o.expr, .name = o.name });
            continue;
        }
        const constraint = if (spill and o.constraint[0] == '=') try std.fmt.allocPrint(map.arena, "=&{s}", .{o.constraint[1..]}) else o.constraint;
        try put(w, "{s}    {s} [{s}] \"{s}\" ({s}),\n", .{ indent, sep, o.name, constraint, try zigExpr(o.expr) });
        sep = " ";
    }
    if (allDropped(outs)) try put(w, "{s}    :\n", .{indent});
    sep = ":";
    var all_ins: std.ArrayList(Operand) = .empty;
    try all_ins.appendSlice(map.arena, ins);
    try all_ins.appendSlice(map.arena, moved.items);
    for (all_ins.items) |o| {
        if (o.dropped) continue;
        // an array decays to its address in C
        const is_local_array = for (locals) |l| {
            if (std.mem.eql(u8, l, o.expr)) break std.mem.indexOf(u8, w.items, try std.fmt.allocPrint(map.arena, "const {s} = [", .{l})) != null;
        } else false;
        const value = if (std.ascii.isDigit(o.expr[0])) try std.fmt.allocPrint(map.arena, "@as(isize, {s})", .{o.expr}) else try zigExpr(o.expr);
        const tie = tieOf(o.constraint);
        const constraint = if (tie) |t| try std.fmt.allocPrint(map.arena, "{d}", .{position[t]}) else o.constraint;
        try put(w, "{s}    {s} [{s}] \"{s}\" ({s}{s}),\n", .{ indent, sep, o.name, constraint, if (is_local_array) "&" else "", value });
        sep = " ";
    }
    if (allDropped(all_ins.items)) try put(w, "{s}    :\n", .{indent});

    // clobbers, at the width of the vector registers the text uses
    try put(w, "{s}    : .{{", .{indent});
    var csep: []const u8 = " ";
    var saw_memory = false;
    var declared_vs: [64]bool = @splat(false);
    var cl = std.mem.tokenizeAny(u8, secs[2], ", \t\r\n");
    while (cl.next()) |tok| {
        const written = std.mem.trim(u8, tok, "\"%");
        const c = if (std.mem.eql(u8, written, "cc")) cc_name else written;
        if (c.len == 0) continue;
        if (std.mem.eql(u8, c, "memory")) saw_memory = true;
        if (std.mem.startsWith(u8, c, "vs")) {
            if (std.fmt.parseInt(usize, c[2..], 10)) |v| declared_vs[v] = true else |_| {}
        }
        if (std.mem.startsWith(u8, c, "xmm")) {
            const n = c[3..];
            const wide = if (std.mem.indexOf(u8, text_out.items, try std.fmt.allocPrint(map.arena, "%%zmm{s}", .{n})) != null) "zmm" else if (std.mem.indexOf(u8, text_out.items, try std.fmt.allocPrint(map.arena, "%%ymm{s}", .{n})) != null) "ymm" else "xmm";
            try put(w, "{s}.{s}{s} = true", .{ csep, wide, n });
        } else try put(w, "{s}.{s} = true", .{ csep, c });
        csep = ", ";
    }
    if (is_power) {
        // the numbered vector registers the text names, as the reference
        // does not declare them
        var seen: [64]bool = @splat(false);
        var lines2 = std.mem.splitScalar(u8, text_out.items, '\n');
        while (lines2.next()) |ln| {
            const t = std.mem.trim(u8, ln, " ");
            if (t.len == 0 or t[0] != 'x' and !std.mem.startsWith(u8, t, "lxv") and !std.mem.startsWith(u8, t, "stxv")) continue;
            var toks = std.mem.tokenizeAny(u8, t, " ,");
            _ = toks.next();
            while (toks.next()) |tok| {
                const v = std.fmt.parseInt(usize, tok, 10) catch continue;
                if (v >= 32 and v < 64) seen[v] = true;
            }
        }
        for (seen, 0..) |on, v| {
            if (!on or declared_vs[v]) continue;
            try put(w, "{s}.vs{d} = true", .{ csep, v });
            csep = ", ";
        }
    }
    // the array-typed memory operands dropped above amount to this
    if (!saw_memory) try put(w, "{s}.memory = true", .{csep});
    try put(w, " }});\n", .{});
}

/// One function: its Zig text.
fn convertFunction(text: []const u8, prec: Precision, conj: bool) !?[]const u8 {
    // the definition, not the prototype
    const def = std.mem.lastIndexOf(u8, text, "static void ") orelse return null;
    const open = std.mem.indexOfScalarPos(u8, text, def, '(').?;
    const close = std.mem.indexOfScalarPos(u8, text, open, ')').?;
    const name = std.mem.trim(u8, text[def + "static void ".len .. open], " ");
    const body_open = std.mem.indexOfScalarPos(u8, text, close, '{') orelse return null;
    const asm_at = std.mem.indexOf(u8, text[body_open..], "__asm__") orelse return null;

    var w: std.ArrayList(u8) = .empty;
    const fname = if (conj) try std.fmt.allocPrint(map.arena, "{c}axpyc{s}", .{ prec.letter, name[1 + "axpy".len ..] }) else name;
    try put(&w, "pub noinline fn {s}(", .{fname});

    var mutable: std.ArrayList([]const u8) = .empty;
    var it = std.mem.splitScalar(u8, text[open + 1 .. close], ',');
    var first = true;
    while (it.next()) |p| {
        const t = std.mem.trim(u8, p, " \t\r\n");
        const star = std.mem.lastIndexOfScalar(u8, t, '*');
        const cut = if (star) |s| s + 1 else std.mem.lastIndexOfScalar(u8, t, ' ').? + 1;
        const pname = std.mem.trim(u8, t[cut..], " ");
        const ctype = std.mem.trim(u8, t[0..cut], " ");
        const ctype_norm = if (star != null) "FLOAT *" else ctype;
        if (!first) try put(&w, ", ", .{});
        first = false;
        // a parameter the assembly writes back is a local copy in Zig
        const written = try isWritten(text[body_open..], pname);
        if (written) try mutable.append(map.arena, pname);
        try put(&w, "{s}{s}: {s}", .{ pname, if (written) "_in" else "", try zigType(ctype_norm, prec.float, pname) });
    }
    try put(&w, ") void {{\n", .{});
    for (mutable.items) |m| try put(&w, "    var {s} = {s}_in;\n", .{ m, m });

    // the statements before the first asm: locals, and possibly the `if`
    // that picks between two blocks
    var locals: std.ArrayList([]const u8) = .empty;
    var condition: ?[]const u8 = null;
    var lines = std.mem.splitScalar(u8, text[body_open + 1 .. body_open + asm_at], '\n');
    while (lines.next()) |raw| {
        const l = std.mem.trim(u8, raw, " \t\r;{");
        if (l.len == 0) continue;
        if (std.mem.startsWith(u8, l, "BLASLONG register ")) {
            const rest = l["BLASLONG register ".len..];
            const eq = std.mem.indexOfScalar(u8, rest, '=').?;
            try put(&w, "    var {s}: isize = {s};\n", .{ std.mem.trim(u8, rest[0..eq], " "), std.mem.trim(u8, rest[eq + 1 ..], " ") });
            try locals.append(map.arena, std.mem.trim(u8, rest[0..eq], " "));
        } else if (std.mem.startsWith(u8, l, "FLOAT ")) {
            const rest = l["FLOAT ".len..];
            const br = std.mem.indexOfScalar(u8, rest, '[').?;
            const rb = std.mem.indexOfScalar(u8, rest, ']').?;
            const vopen = std.mem.indexOfScalar(u8, rest, '{').?;
            const vclose = std.mem.lastIndexOfScalar(u8, rest, '}').?;
            try put(&w, "    const {s} = [{s}]{s}{{{s}}};\n", .{ rest[0..br], rest[br + 1 .. rb], prec.float, rest[vopen + 1 .. vclose] });
            try locals.append(map.arena, rest[0..br]);
        } else if (std.mem.startsWith(u8, l, "static const ")) {
            // `static const double mvec[2] = { -1.0, 1.0 };`
            const rest = l["static const ".len..];
            const sp = std.mem.indexOfScalar(u8, rest, ' ').?;
            const br = std.mem.indexOfScalar(u8, rest, '[').?;
            const rb = std.mem.indexOfScalar(u8, rest, ']').?;
            const vopen = std.mem.indexOfScalar(u8, rest, '{').?;
            const vclose = std.mem.lastIndexOfScalar(u8, rest, '}').?;
            try put(&w, "    const {s} = [{s}]{s}{{{s}}};\n", .{ rest[sp + 1 .. br], rest[br + 1 .. rb], prec.float, rest[vopen + 1 .. vclose] });
            try locals.append(map.arena, rest[sp + 1 .. br]);
        } else if (std.mem.startsWith(u8, l, "const double *") or std.mem.startsWith(u8, l, "const float *")) {
            const rest = l[std.mem.indexOfScalar(u8, l, '*').? + 1 ..];
            const eq = std.mem.indexOfScalar(u8, rest, '=').?;
            try put(&w, "    const {s}: [*]const {s} = &{s};\n", .{ std.mem.trim(u8, rest[0..eq], " "), prec.float, std.mem.trim(u8, rest[eq + 1 ..], " ") });
        } else if (std.mem.startsWith(u8, l, "__vector unsigned char ")) {
            const rest = l["__vector unsigned char ".len..];
            const eq = std.mem.indexOfScalar(u8, rest, '=').?;
            try put(&w, "    const {s}: @Vector(16, u8) = .{s};\n", .{ std.mem.trim(u8, rest[0..eq], " "), std.mem.trim(u8, rest[eq + 1 ..], " ") });
        } else if (std.mem.startsWith(u8, l, "__vector ")) {
            const rest = l["__vector ".len..];
            const sp = std.mem.indexOfScalar(u8, rest, ' ').?;
            const lanes: usize = if (std.mem.eql(u8, rest[0..sp], "double")) 2 else 4;
            const decl = std.mem.trim(u8, rest[sp + 1 ..], " ");
            if (std.mem.indexOfScalar(u8, decl, '=')) |eq| {
                try put(&w, "    const {s}: @Vector({d}, {s}) = .{s};\n", .{ std.mem.trim(u8, decl[0..eq], " "), lanes, prec.float, std.mem.trim(u8, decl[eq + 1 ..], " ") });
            } else try put(&w, "    var {s}: @Vector({d}, {s}) = undefined;\n", .{ decl, lanes, prec.float });
        } else if (std.mem.startsWith(u8, l, "long ")) {
            try put(&w, "    var {s}: isize = undefined;\n", .{std.mem.trim(u8, l["long ".len..], " ")});
        } else if (std.mem.startsWith(u8, l, "if")) {
            const po = std.mem.indexOfScalar(u8, l, '(').?;
            const pc = std.mem.lastIndexOfScalar(u8, l, ')').?;
            condition = std.mem.trim(u8, l[po + 1 .. pc], " ");
        } else map.die("cannot convert the statement `{s}` in {s}", .{ l, name });
    }

    // every asm statement, by paren matching
    var stmts: std.ArrayList([]const u8) = .empty;
    var pos = body_open;
    while (std.mem.indexOfPos(u8, text, pos, "__asm__")) |at| {
        const po = std.mem.indexOfScalarPos(u8, text, at, '(').?;
        var depth: usize = 0;
        var in_str = false;
        var e = po;
        while (e < text.len) : (e += 1) {
            const c = text[e];
            if (c == '"' and text[e - 1] != '\\') in_str = !in_str;
            if (in_str) continue;
            if (c == '(') depth += 1;
            if (c == ')') {
                depth -= 1;
                if (depth == 0) break;
            }
        }
        try stmts.append(map.arena, text[po + 1 .. e]);
        pos = e;
    }
    if (stmts.items.len == 2 and condition == null) map.die("two asm blocks with no condition in {s}", .{name});

    for (stmts.items, 0..) |stmt, k| {
        const indent: []const u8 = if (stmts.items.len == 2) "        " else "    ";
        if (stmts.items.len == 2) {
            if (k == 0) try put(&w, "    if ({s}) {{\n", .{condition.?}) else try put(&w, "    }} else {{\n", .{});
        }
        try emitAsm(&w, indent, stmt, locals.items, name);
    }
    if (stmts.items.len == 2) try put(&w, "    }}\n", .{});
    try put(&w, "}}\n", .{});
    return w.items;
}

/// Converts one microkernel file into `out_dir`.
/// An architecture's block kernels: where they are in OpenBLAS, and where
/// they go here.
const Target = struct {
    /// Under `kernel/`.
    ob_dir: []const u8,
    /// Under `blas1/axpy/`.
    out: []const u8,
    /// Whether the block is a function in the driver's own file, not in an
    /// included `_microk_` header chosen by an `#if` chain.
    in_file: bool,
    /// Whether it is POWER, which is preprocessed for a VSX target.
    power: bool = false,
    /// Whether there is a driver `#if` chain to index.
    index: bool = true,
};

const targets = [_]Target{
    .{ .ob_dir = "x86_64", .out = "x86_64", .in_file = false },
    .{ .ob_dir = "zarch", .out = "s390x", .in_file = true },
    .{ .ob_dir = "power", .out = "power", .in_file = false, .power = true, .index = false },
};

/// Converts one file of blocks into `out_dir`.
fn convertFile(t: Target, path: []const u8, out_dir: []const u8) !void {
    const base = std.fs.path.basename(path);
    const prec = precisionOf(base[0]);
    const stem = base[0 .. base.len - 2];
    const file = if (t.in_file)
        try std.fmt.allocPrint(map.arena, "{s}.zig", .{stem})
    else blk: {
        // `daxpy_microk_haswell-2` is `daxpy_haswell`
        const core = stem[std.mem.indexOf(u8, stem, "_microk_").? + "_microk_".len ..];
        const core_name = if (std.mem.endsWith(u8, core, "-2")) core[0 .. core.len - 2] else core;
        break :blk try std.fmt.allocPrint(map.arena, "{c}axpy_{s}.zig", .{ prec.letter, core_name });
    };

    const original = map.read(path).?;
    const has_conj = prec.complex and std.mem.indexOf(u8, original, "CONJ") != null;
    var out: std.ArrayList(u8) = .empty;
    try put(&out,
        \\//! `{c}axpy` block kernel from OpenBLAS `kernel/{s}/{s}`, the assembly
        \\//! as written. Generated by `tools/blas_asm/microk.zig`; re-derive it
        \\//! rather than editing it.
        \\
        \\
    , .{ prec.letter, t.ob_dir, base });
    var wrote = false;
    for ([_]bool{ false, true }) |conj| {
        if (conj and !has_conj) continue;
        const text = try runCpp(path, conj);
        if (try convertFunction(text, prec, conj)) |z| {
            try put(&out, "{s}\n", .{z});
            wrote = true;
        }
    }
    if (!wrote) {
        std.debug.print("skipped {s}: no assembly\n", .{base});
        return;
    }
    const dest = try map.path(&.{ out_dir, file });
    try Io.Dir.cwd().writeFile(map.io, .{ .sub_path = dest, .data = out.items });
    _ = try std.process.run(map.arena, map.io, .{ .argv = &.{ "zig", "fmt", dest } });
    std.debug.print("wrote {s}\n", .{dest});
}

pub fn main(init: std.process.Init) !void {
    const a = init.arena.allocator();
    const checkout = init.environ_map.get("OPENBLAS") orelse
        map.die("set OPENBLAS to the OpenBLAS checkout", .{});
    try map.init(a, init.io, checkout);
    if (init.environ_map.get("BLAS_OUT")) |o| map.out_root = o;

    for (targets) |t| {
        cc_name = if (std.mem.eql(u8, t.out, "s390x")) "ps" else "cc";
        is_power = t.power;
        cpp_extra = if (t.power) &.{
            "-target",                                                           "powerpc64le-linux",
            "-mcpu=pwr10",                                                       "-DPOWER10",
            "-D__VEC__",                                                         "-DXXSPLTD_S(T,A,z)=\"xxspltd \" #T \", \" #A \", \" #z \"\\n\\t\"",
            "-DXXSWAPD_S(T,A)=\"xxswapd \" #T \", \" #A \"\\n\\t\"",             "-DXXMRGHD_S(T,A,B)=\"xxmrghd \" #T \", \" #A \", \" #B \"\\n\\t\"",
            "-DXXMRGLD_S(T,A,B)=\"xxmrgld \" #T \", \" #A \", \" #B \"\\n\\t\"",
        } else &.{};
        const out_dir = try std.fmt.allocPrint(a, "{s}/{s}/microk", .{ map.out_root, t.out });
        try Io.Dir.cwd().createDirPath(init.io, out_dir);

        var dir = try Io.Dir.cwd().openDir(init.io, try map.path(&.{ map.root, "kernel", t.ob_dir }), .{ .iterate = true });
        defer dir.close(init.io);
        var names: std.ArrayList([]const u8) = .empty;
        var it = dir.iterate();
        while (try it.next(init.io)) |e| {
            if (e.kind != .file or !std.mem.endsWith(u8, e.name, ".c")) continue;
            const wanted = if (t.in_file)
                std.mem.eql(u8, e.name[1..], "axpy.c") and std.mem.indexOfScalar(u8, "sdcz", e.name[0]) != null
            else
                std.mem.indexOf(u8, e.name, "axpy_microk_") != null;
            if (wanted) try names.append(a, try a.dupe(u8, e.name));
        }
        std.mem.sort([]const u8, names.items, {}, struct {
            fn lt(_: void, x: []const u8, y: []const u8) bool {
                return std.mem.lessThan(u8, x, y);
            }
        }.lt);
        for (names.items) |n| try convertFile(t, try map.path(&.{ map.root, "kernel", t.ob_dir, n }), out_dir);
        if (t.index) try writeIndex(a, t, out_dir);
    }
}

/// One `#if`/`#elif` of a driver's chain: the cores it names, and the block
/// it includes.
const Link = struct { cores: []const []const u8, block: []const u8 };

/// The chain at the top of a driver, in order.
fn readChain(a: std.mem.Allocator, text: []const u8) ![]const Link {
    var out: std.ArrayList(Link) = .empty;
    var cores: []const []const u8 = &.{};
    var it = std.mem.splitScalar(u8, text, '\n');
    while (it.next()) |raw| {
        const line = std.mem.trim(u8, raw, " \t\r");
        const is_if = std.mem.startsWith(u8, line, "#if defined") or std.mem.startsWith(u8, line, "#elif defined");
        if (is_if) {
            var list: std.ArrayList([]const u8) = .empty;
            const cond = line[std.mem.indexOf(u8, line, "defined").?..];
            var parts = std.mem.splitSequence(u8, cond, "||");
            while (parts.next()) |p| {
                const t = std.mem.trim(u8, p, " ");
                const name = std.mem.trim(u8, t["defined".len..], " ()");
                try list.append(a, name);
            }
            cores = list.items;
        } else if (std.mem.startsWith(u8, line, "#include") and std.mem.indexOf(u8, line, "_microk_") != null) {
            const q1 = std.mem.indexOfScalar(u8, line, '"').?;
            const q2 = std.mem.lastIndexOfScalar(u8, line, '"').?;
            const file = line[q1 + 1 .. q2];
            const core = file[std.mem.indexOf(u8, file, "_microk_").? + "_microk_".len ..];
            const stem = if (std.mem.endsWith(u8, core, "-2.c")) core[0 .. core.len - 4] else core;
            try out.append(a, .{ .cores = cores, .block = stem });
        } else if (std.mem.startsWith(u8, line, "#ifndef HAVE_KERNEL")) break;
    }
    return out.items;
}

/// The block size a driver hands its block: the first `n & -N` it computes.
fn readMask(text: []const u8) []const u8 {
    const key = "n & -";
    const at = std.mem.indexOf(u8, text, key).? + key.len;
    var e = at;
    while (e < text.len and std.ascii.isDigit(text[e])) e += 1;
    return text[at..e];
}

/// The kernels of a converted block, by role: the first function is the real
/// or non-conjugated one, a second the conjugated.
fn functionsOf(a: std.mem.Allocator, path: []const u8) ![]const []const u8 {
    const text = map.read(path) orelse return &.{};
    var out: std.ArrayList([]const u8) = .empty;
    var it = std.mem.splitScalar(u8, text, '\n');
    while (it.next()) |line| {
        const key = "pub noinline fn ";
        if (!std.mem.startsWith(u8, line, key)) continue;
        const rest = line[key.len..];
        try out.append(a, rest[0..std.mem.indexOfScalar(u8, rest, '(').?]);
    }
    return out.items;
}

/// The kernels a driver serves: real, or complex and conjugated.
fn rolesOf(a: std.mem.Allocator, letter: u8) ![]const []const u8 {
    const out = try a.alloc([]const u8, if (letter == 'c' or letter == 'z') 2 else 1);
    out[0] = try std.fmt.allocPrint(a, "{c}axpy", .{letter});
    if (out.len == 2) out[1] = try std.fmt.allocPrint(a, "{c}axpyc", .{letter});
    return out;
}

/// Which block each core's driver uses, and at what block size.
fn writeIndex(a: std.mem.Allocator, t: Target, out_dir: []const u8) !void {
    const arch = t.ob_dir;
    var w: std.ArrayList(u8) = .empty;
    try put(&w,
        \\//! Which block kernel each core's C driver uses, read from the `#if`
        \\//! chain at the top of `kernel/x86_64/<p>axpy.c`. Generated by
        \\//! `tools/blas_asm/microk.zig`; re-derive it rather than editing it.
        \\
        \\const std = @import("std");
        \\
        \\/// What the driver calls for the contiguous part.
        \\pub const Block = enum {{
        \\    /// An assembly function of `module`.
        \\    assembly,
        \\    /// C intrinsics upstream, written on vectors in `kernel.zig`.
        \\    vector,
        \\    /// Nothing: the core's driver has no block, the scalar loops do it all.
        \\    none,
        \\}};
        \\
        \\pub const Entry = struct {{
        \\    role: []const u8,
        \\    cores: []const []const u8,
        \\    /// The driver's block size, `n & -mask`.
        \\    mask: comptime_int,
        \\    block: Block,
        \\    /// The block's function, when it is assembly.
        \\    kernel: []const u8,
        \\    module: type,
        \\}};
        \\
        \\pub const entries = [_]Entry{{
        \\
    , .{});

    const all = (try map.cores(arch)).?;
    for ([_]u8{ 's', 'd', 'c', 'z' }) |letter| {
        const driver = try std.fmt.allocPrint(a, "{c}axpy.c", .{letter});
        const text = map.read(try map.path(&.{ map.root, "kernel", arch, driver })).?;
        // the block is either the one function in the driver's own file, for
        // every core that uses the driver, or what the `#if` chain includes
        const chain: []const Link = if (t.in_file) &.{.{ .cores = all, .block = "" }} else try readChain(a, text);
        const mask = readMask(text);
        const upper = try std.fmt.allocPrint(a, "{c}", .{std.ascii.toUpper(letter)});
        var named_cores: std.ArrayList([]const u8) = .empty;
        for (chain) |link| {
            // the cores whose kernel list selects this driver and whose macro the link names
            var cores: std.ArrayList(u8) = .empty;
            for (all) |core| {
                const uses = for (try map.resolve(arch, core, upper)) |s| {
                    if (std.mem.eql(u8, s, driver)) break true;
                } else false;
                if (!uses) continue;
                const named = for (link.cores) |c| {
                    if (std.mem.eql(u8, c, core)) break true;
                } else false;
                if (!named) continue;
                const lower = try std.ascii.allocLowerString(a, core);
                try cores.print(a, " \"{s}\",", .{lower});
                try named_cores.append(a, core);
            }
            if (cores.items.len == 0) continue;

            const block_file = if (t.in_file)
                try std.fmt.allocPrint(a, "{c}axpy.zig", .{letter})
            else
                try std.fmt.allocPrint(a, "{c}axpy_{s}.zig", .{ letter, link.block });
            const fns = try functionsOf(a, try map.path(&.{ out_dir, block_file }));
            for (try rolesOf(a, letter), 0..) |role, k| {
                const kernel = if (fns.len == 0) "" else fns[@min(k, fns.len - 1)];
                const module = if (fns.len == 0) "void" else try std.fmt.allocPrint(a, "@import(\"{s}\")", .{block_file});
                try put(&w, "    .{{ .role = \"{s}\", .cores = &.{{{s} }}, .mask = {s}, .block = .{s}, .kernel = \"{s}\", .module = {s} }},\n", .{ role, cores.items, mask, if (fns.len == 0) "vector" else "assembly", kernel, module });
            }
        }
        // the cores that use the driver but no link names: its own scalar fallback
        var rest: std.ArrayList(u8) = .empty;
        for (all) |core| {
            const uses = for (try map.resolve(arch, core, upper)) |s| {
                if (std.mem.eql(u8, s, driver)) break true;
            } else false;
            if (!uses) continue;
            const named = for (named_cores.items) |c| {
                if (std.mem.eql(u8, c, core)) break true;
            } else false;
            if (named) continue;
            try rest.print(a, " \"{s}\",", .{try std.ascii.allocLowerString(a, core)});
        }
        if (rest.items.len > 0) for (try rolesOf(a, letter)) |role| {
            try put(&w, "    .{{ .role = \"{s}\", .cores = &.{{{s} }}, .mask = {s}, .block = .none, .kernel = \"\", .module = void }},\n", .{ role, rest.items, mask });
        };
    }
    try put(&w, "}};\n\n", .{});
    try w.appendSlice(a,
        \\/// The entry for `role` on `core`, if its driver is C.
        \\pub fn find(comptime role: []const u8, comptime core: []const u8) ?Entry {
        \\    @setEvalBranchQuota(1_000_000);
        \\    inline for (entries) |e| {
        \\        if (comptime std.mem.eql(u8, e.role, role)) {
        \\            inline for (e.cores) |c| {
        \\                if (comptime std.mem.eql(u8, c, core)) return e;
        \\            }
        \\        }
        \\    }
        \\    return null;
        \\}
        \\
    );
    const dest = try map.path(&.{ out_dir, "__root__.zig" });
    try Io.Dir.cwd().writeFile(map.io, .{ .sub_path = dest, .data = w.items });
    _ = try std.process.run(a, map.io, .{ .argv = &.{ "zig", "fmt", dest } });
    std.debug.print("wrote {s}\n", .{dest});
}
