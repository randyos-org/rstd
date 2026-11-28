//! TrueType instruction opcodes, their lengths and stack effects.
//!
//! Names from fontations `read-fonts/src/tables/glyf/bytecode/opcode.rs`
//! (commit 4f115482, MIT); the pop and push counts are the OpenType
//! specification's, as FreeType 2.14.3 `ttinterp.c` `Pop_Push_Count` has them.

/// An opcode; the PUSHB, PUSHW, MDRP and MIRP families and the undefined
/// opcodes have no names (`isPush`, `mdrp_first`, `mirp_first`).
pub const Opcode = enum(u8) {
    svtca0 = 0x00,
    svtca1 = 0x01,
    spvtca0 = 0x02,
    spvtca1 = 0x03,
    sfvtca0 = 0x04,
    sfvtca1 = 0x05,
    spvtl0 = 0x06,
    spvtl1 = 0x07,
    sfvtl0 = 0x08,
    sfvtl1 = 0x09,
    spvfs = 0x0A,
    sfvfs = 0x0B,
    gpv = 0x0C,
    gfv = 0x0D,
    sfvtpv = 0x0E,
    isect = 0x0F,
    srp0 = 0x10,
    srp1 = 0x11,
    srp2 = 0x12,
    szp0 = 0x13,
    szp1 = 0x14,
    szp2 = 0x15,
    szps = 0x16,
    sloop = 0x17,
    rtg = 0x18,
    rthg = 0x19,
    smd = 0x1A,
    @"else" = 0x1B,
    jmpr = 0x1C,
    scvtci = 0x1D,
    sswci = 0x1E,
    ssw = 0x1F,
    dup = 0x20,
    pop = 0x21,
    clear = 0x22,
    swap = 0x23,
    depth = 0x24,
    cindex = 0x25,
    mindex = 0x26,
    alignpts = 0x27,
    utp = 0x29,
    loopcall = 0x2A,
    call = 0x2B,
    fdef = 0x2C,
    endf = 0x2D,
    mdap0 = 0x2E,
    mdap1 = 0x2F,
    iup0 = 0x30,
    iup1 = 0x31,
    shp0 = 0x32,
    shp1 = 0x33,
    shc0 = 0x34,
    shc1 = 0x35,
    shz0 = 0x36,
    shz1 = 0x37,
    shpix = 0x38,
    ip = 0x39,
    msirp0 = 0x3A,
    msirp1 = 0x3B,
    alignrp = 0x3C,
    rtdg = 0x3D,
    miap0 = 0x3E,
    miap1 = 0x3F,
    npushb = 0x40,
    npushw = 0x41,
    ws = 0x42,
    rs = 0x43,
    wcvtp = 0x44,
    rcvt = 0x45,
    gc0 = 0x46,
    gc1 = 0x47,
    scfs = 0x48,
    md0 = 0x49,
    md1 = 0x4A,
    mppem = 0x4B,
    mps = 0x4C,
    flipon = 0x4D,
    flipoff = 0x4E,
    debug = 0x4F,
    lt = 0x50,
    lteq = 0x51,
    gt = 0x52,
    gteq = 0x53,
    eq = 0x54,
    neq = 0x55,
    odd = 0x56,
    even = 0x57,
    @"if" = 0x58,
    eif = 0x59,
    @"and" = 0x5A,
    @"or" = 0x5B,
    not = 0x5C,
    deltap1 = 0x5D,
    sdb = 0x5E,
    sds = 0x5F,
    add = 0x60,
    sub = 0x61,
    div = 0x62,
    mul = 0x63,
    abs = 0x64,
    neg = 0x65,
    floor = 0x66,
    ceiling = 0x67,
    round00 = 0x68,
    round01 = 0x69,
    round10 = 0x6A,
    round11 = 0x6B,
    nround00 = 0x6C,
    nround01 = 0x6D,
    nround10 = 0x6E,
    nround11 = 0x6F,
    wcvtf = 0x70,
    deltap2 = 0x71,
    deltap3 = 0x72,
    deltac1 = 0x73,
    deltac2 = 0x74,
    deltac3 = 0x75,
    sround = 0x76,
    s45round = 0x77,
    jrot = 0x78,
    jrof = 0x79,
    roff = 0x7A,
    rutg = 0x7C,
    rdtg = 0x7D,
    sangw = 0x7E,
    aa = 0x7F,
    flippt = 0x80,
    fliprgon = 0x81,
    fliprgoff = 0x82,
    scanctrl = 0x85,
    sdpvtl0 = 0x86,
    sdpvtl1 = 0x87,
    getinfo = 0x88,
    idef = 0x89,
    roll = 0x8A,
    max = 0x8B,
    min = 0x8C,
    scantype = 0x8D,
    instctrl = 0x8E,
    getvariation = 0x91,
    getdata = 0x92,
    _,
};

pub const pushb_first = 0xB0;
pub const pushw_first = 0xB8;
pub const mdrp_first = 0xC0;
pub const mirp_first = 0xE0;

/// NPUSHB, NPUSHW, PUSHB[n] and PUSHW[n].
pub fn isPush(op: u8) bool {
    return op == 0x40 or op == 0x41 or (op >= pushb_first and op < mdrp_first);
}

/// NPUSHW and PUSHW[n]: 16-bit operands.
pub fn isPushWords(op: u8) bool {
    return op == 0x41 or (op >= pushw_first and op < mdrp_first);
}

/// Each opcode's length in bytes; negative for NPUSHB (-1) and NPUSHW (-2):
/// two bytes plus the count byte's operands of that size.
pub const lengths: [256]i8 = blk: {
    var t: [256]i8 = @splat(1);
    t[0x40] = -1;
    t[0x41] = -2;
    for (0..8) |n| {
        t[pushb_first + n] = 2 + n;
        t[pushw_first + n] = 3 + 2 * n;
    }
    break :blk t;
};

/// How many values an instruction pops and pushes; the variable ones (the
/// pushes, GETVARIATION, the loops and the DELTAs' pairs) count them itself.
pub const StackEffect = struct { pop: u8, push: u8 };

pub const stack_effects: [256]StackEffect = blk: {
    var t: [256]StackEffect = @splat(.{ .pop = 0, .push = 0 });
    const set = struct {
        fn f(table: *[256]StackEffect, ops: []const u8, pop: u8, push: u8) void {
            for (ops) |op| table[op] = .{ .pop = pop, .push = push };
        }
    }.f;
    // one popped
    set(&t, &.{ 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x1A, 0x1C, 0x1D, 0x1E, 0x1F, 0x21, 0x26, 0x29, 0x2B, 0x2C, 0x2E, 0x2F, 0x34, 0x35, 0x36, 0x37, 0x38, 0x4F, 0x58, 0x5D, 0x5E, 0x5F, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77, 0x7E, 0x7F, 0x85, 0x89, 0x8D }, 1, 0);
    // two popped
    set(&t, &.{ 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0x27, 0x2A, 0x3A, 0x3B, 0x3E, 0x3F, 0x42, 0x44, 0x48, 0x70, 0x78, 0x79, 0x81, 0x82, 0x86, 0x87, 0x8E }, 2, 0);
    // one in, one out
    set(&t, &.{ 0x25, 0x43, 0x45, 0x46, 0x47, 0x56, 0x57, 0x5C, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69, 0x6A, 0x6B, 0x6C, 0x6D, 0x6E, 0x6F, 0x88 }, 1, 1);
    // two in, one out
    set(&t, &.{ 0x49, 0x4A, 0x50, 0x51, 0x52, 0x53, 0x54, 0x55, 0x5A, 0x5B, 0x60, 0x61, 0x62, 0x63, 0x8B, 0x8C }, 2, 1);
    // pushes only
    set(&t, &.{ 0x24, 0x4B, 0x4C, 0x92 }, 0, 1);
    set(&t, &.{ 0x0C, 0x0D }, 0, 2);
    t[0x0F] = .{ .pop = 5, .push = 0 }; // ISECT
    t[0x20] = .{ .pop = 1, .push = 2 }; // DUP
    t[0x23] = .{ .pop = 2, .push = 2 }; // SWAP
    t[0x8A] = .{ .pop = 3, .push = 3 }; // ROLL
    for (0..8) |n| {
        t[pushb_first + n] = .{ .pop = 0, .push = n + 1 };
        t[pushw_first + n] = .{ .pop = 0, .push = n + 1 };
    }
    for (mdrp_first..mirp_first) |op| t[op] = .{ .pop = 1, .push = 0 };
    for (mirp_first..256) |op| t[op] = .{ .pop = 2, .push = 0 };
    break :blk t;
};
