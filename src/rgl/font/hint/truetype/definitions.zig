//! Function and instruction definitions: FreeType 2.14.3 `ttinterp.c`'s
//! `FDefs` and `IDefs` tables, records in definition order found by
//! number, with `maxFunc` and `maxIns`.

const Program = @import("program.zig").Program;

pub const Definition = struct {
    program: Program = .font,
    /// The offset after FDEF or IDEF.
    start: usize = 0,
    /// The offset of its ENDF.
    end: usize = 0,
    /// The function number, or the opcode.
    number: u32 = 0,
    active: bool = false,
};

pub const Table = struct {
    /// Room for `maxp`'s count.
    defs: []Definition,
    count: usize = 0,
    /// The largest number defined (`maxFunc`, `maxIns`).
    max_number: u32 = 0,

    pub fn reset(t: *Table) void {
        t.count = 0;
        t.max_number = 0;
    }

    /// The definition numbered `n`, defined or not.
    pub fn find(t: *Table, n: u64) ?*Definition {
        for (t.defs[0..t.count]) |*d| if (d.number == n) return d;
        return null;
    }

    /// The record for `n`: its old one, else a new one, null past room.
    pub fn slot(t: *Table, n: u64) ?*Definition {
        if (t.find(n)) |d| return d;
        if (t.count >= t.defs.len) return null;
        t.count += 1;
        return &t.defs[t.count - 1];
    }
};
