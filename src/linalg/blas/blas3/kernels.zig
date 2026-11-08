//! The level-3 kernel registry: the modules with kernels over the gemm tile
//! (`sets`), their portable install, and each architecture's per-core
//! kernels (`Registry`, for `kernels/<arch>.zig`).

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../types.zig");
const util = @import("../util.zig");
const Vtable = @import("../Vtable.zig");
const gk = @import("gemm/kernel.zig");

const Slot = util.Slot;
const XDouble = types.XDouble;
const isComplex = types.isComplex;

/// Every module with kernels over the gemm tile: each `gk.Kernels(T, geo)` has
/// a public function per vtable slot it fills, named as the slot less its
/// precision letter.
pub const sets = .{ gk, @import("symm/kernel.zig"), @import("trmm/kernel.zig"), @import("trsm/kernel.zig"), @import("../ext/neg/kernel.zig"), @import("../ext/laswp/kernel.zig"), three_m, @import("gemm/small_kernel.zig") };

/// The 3M kernels, whose tile is the real one of the complex type's parts:
/// their unrolls are their own (`<p>gemm3m_unroll_m`, `_n`).
const three_m = @import("gemm3m/kernel.zig");

/// The precisions these kernels serve, by vtable prefix.
const precisions = .{ .{ "s", f32 }, .{ "d", f64 }, .{ "q", XDouble }, .{ "c", types.Complex32 }, .{ "z", types.Complex64 }, .{ "x", types.ComplexXD } };
/// Those a core's object may replace: `q` and `x` slots keep the portable kernels.
const core_precisions = .{ precisions[0], precisions[1], precisions[3], precisions[4] };

/// The slot that says whose kernels a precision has.
fn witness(comptime T: type) []const u8 {
    return if (isComplex(T)) "gemm_kernel_n" else "gemm_kernel";
}

/// The names of `K`'s functions that are slots with prefix `pre`.
fn slotNames(comptime pre: []const u8, comptime K: type) []const []const u8 {
    comptime {
        var out: []const []const u8 = &.{};
        for (std.meta.declarations(K)) |d| {
            if (@typeInfo(@TypeOf(@field(K, d))) == .@"fn" and @hasField(Vtable, pre ++ d)) out = out ++ [_][]const u8{d};
        }
        return out;
    }
}

/// The portable kernel of `T` named `name` (a slot less its letter), from
/// whichever of `sets` has it: `portable.kernel`'s source for level 3.
pub fn find(comptime T: type, comptime name: []const u8) @TypeOf(&@field(Owner(T, name).Kernels(T, gk.geometry(T)), name)) {
    return &@field(Owner(T, name).Kernels(T, gk.geometry(T)), name);
}

fn Owner(comptime T: type, comptime name: []const u8) type {
    inline for (sets) |s| {
        if (@hasDecl(s.Kernels(T, gk.geometry(T)), name)) return s;
    }
    @compileError("no level-3 kernel " ++ name ++ " for " ++ @typeName(T));
}

/// Points `vt`'s level-3 slots at the kernels of each of `list`, and its
/// gemm unrolls at the tiles, through `source` (`function`, `unroll`).
fn fill(vt: *Vtable, comptime list: anytype, comptime source: anytype) void {
    inline for (list) |p| {
        inline for (sets) |s| {
            // `q` and `x` alias `d` and `z` where `XDouble` is `f64`, and
            // nothing reaches their slots but gemm's (`portable.letter`)
            const alias = XDouble == f64 and (p[0][0] == 'q' or p[0][0] == 'x');
            if (comptime alias and s != gk) continue;
            inline for (comptime slotNames(p[0], s.Kernels(p[1], gk.geometry(p[1])))) |name| {
                @field(vt, p[0] ++ name) = source.function(p[0], p[1], s, name);
            }
        }
        const um = source.unroll(p[0], p[1], "m");
        const un = source.unroll(p[0], p[1], "n");
        @field(vt, p[0] ++ "gemm_unroll_m") = um;
        @field(vt, p[0] ++ "gemm_unroll_n") = un;
        @field(vt, p[0] ++ "gemm_unroll_mn") = @max(um, un);
        if (comptime isComplex(p[1]) and !(XDouble == f64 and p[0][0] == 'x')) {
            setUnrolls3m(vt, p[0], source.unroll3m(p[0], p[1], "m"), source.unroll3m(p[0], p[1], "n"));
        }
    }
}

fn setUnrolls3m(vt: *Vtable, comptime pre: []const u8, um: i32, un: i32) void {
    @field(vt, pre ++ "gemm3m_unroll_m") = um;
    @field(vt, pre ++ "gemm3m_unroll_n") = un;
    @field(vt, pre ++ "gemm3m_unroll_mn") = @max(um, un);
}

/// Puts back the portable kernels' unrolls where `vt`'s slots still hold
/// them: a core's blocking (`params.installCore`), set after them, writes
/// `param.h`'s, and no per-core object may come to replace both.
pub fn settle(vt: *Vtable) void {
    inline for (precisions) |p| {
        const K = gk.Kernels(p[1], gk.geometry(p[1]));
        if (@field(vt, p[0] ++ witness(p[1])) == &@field(K, witness(p[1]))) {
            @field(vt, p[0] ++ "gemm_unroll_m") = K.unroll_m;
            @field(vt, p[0] ++ "gemm_unroll_n") = K.unroll_n;
            @field(vt, p[0] ++ "gemm_unroll_mn") = @max(K.unroll_m, K.unroll_n);
        }
        if (comptime isComplex(p[1])) {
            const K3 = three_m.Kernels(p[1], gk.geometry(p[1]));
            if (@field(vt, p[0] ++ "gemm3m_kernel") == &K3.gemm3m_kernel) setUnrolls3m(vt, p[0], K3.unroll_m, K3.unroll_n);
        }
    }
}

/// The portable kernels: the main library's own build of these.
pub fn installGeneric(vt: *Vtable) void {
    fill(vt, precisions, struct {
        fn function(comptime pre: []const u8, comptime T: type, comptime s: type, comptime name: []const u8) @typeInfo(@FieldType(Vtable, pre ++ name)).optional.child {
            return &@field(s.Kernels(T, gk.geometry(T)), name);
        }
        fn unroll(comptime _: []const u8, comptime T: type, comptime which: []const u8) i32 {
            return @field(gk.Kernels(T, gk.geometry(T)), "unroll_" ++ which);
        }
        fn unroll3m(comptime _: []const u8, comptime T: type, comptime which: []const u8) i32 {
            return @field(three_m.Kernels(T, gk.geometry(T)), "unroll_" ++ which);
        }
    });
}

/// Points `vt`'s slots for `T` (prefix `pre`) at the kernels of tile `geo`,
/// and its unrolls at the tile: for tests, which run the drivers over tiles
/// this machine does not choose.
pub fn installTile(vt: *Vtable, comptime pre: []const u8, comptime T: type, comptime geo: gk.Geometry) void {
    fill(vt, .{.{ pre, T }}, struct {
        fn function(comptime p: []const u8, comptime U: type, comptime s: type, comptime name: []const u8) @typeInfo(@FieldType(Vtable, p ++ name)).optional.child {
            return &@field(s.Kernels(U, geo), name);
        }
        fn unroll(comptime _: []const u8, comptime U: type, comptime which: []const u8) i32 {
            return @field(gk.Kernels(U, geo), "unroll_" ++ which);
        }
        fn unroll3m(comptime _: []const u8, comptime U: type, comptime which: []const u8) i32 {
            return @field(three_m.Kernels(U, geo), "unroll_" ++ which);
        }
    });
}

/// Whether a core's features can change this target's tiles (`vectorBits`):
/// where none can (SPARC), a core's object would only hold the portable
/// kernels again, so it holds none of these.
fn tilesVary() bool {
    return switch (sysinfo.target.cpu.arch) {
        .sparc, .sparc64 => false,
        else => true,
    };
}

/// The per-core kernels, for an architecture's registry (`kernels/<arch>.zig`):
/// every per-core object exports them, built for its own features, with
/// their unrolls as data (the tile depends on the object's features, which
/// the main library does not have). They never make a core need an object.
pub fn Registry(comptime Core: type) type {
    return struct {
        pub fn has(comptime core: Core) bool {
            _ = core;
            return false;
        }

        fn Unrolls(comptime T: type) type {
            return struct {
                const m: i32 = gk.Kernels(T, gk.geometry(T)).unroll_m;
                const n: i32 = gk.Kernels(T, gk.geometry(T)).unroll_n;
            };
        }

        fn Unrolls3m(comptime T: type) type {
            return struct {
                const m: i32 = three_m.Kernels(T, gk.geometry(T)).unroll_m;
                const n: i32 = three_m.Kernels(T, gk.geometry(T)).unroll_n;
            };
        }

        pub fn exportCore(comptime core: Core) void {
            if (comptime !tilesVary()) return;
            inline for (core_precisions) |p| {
                inline for (sets) |s| {
                    // a set may hand a core upstream's own kernels
                    const S = if (@hasDecl(s, "forCore")) s.forCore(core) else s;
                    const K = S.Kernels(p[1], gk.geometry(p[1]));
                    inline for (comptime slotNames(p[0], K)) |name| {
                        @export(&@field(K, name), .{ .name = util.kernelSymbol(@field(Slot, p[0] ++ name), core), .visibility = .hidden });
                    }
                }
                inline for (.{ "m", "n" }) |w| {
                    @export(&@field(Unrolls(p[1]), w), .{ .name = util.kernelSymbol(@field(Slot, p[0] ++ "gemm_unroll_" ++ w), core), .visibility = .hidden });
                    if (comptime isComplex(p[1])) @export(&@field(Unrolls3m(p[1]), w), .{ .name = util.kernelSymbol(@field(Slot, p[0] ++ "gemm3m_unroll_" ++ w), core), .visibility = .hidden });
                }
            }
        }

        pub fn install(vt: *Vtable, comptime core: Core) void {
            if (comptime !tilesVary()) return;
            fill(vt, core_precisions, struct {
                fn function(comptime pre: []const u8, comptime _: type, comptime _: type, comptime name: []const u8) @typeInfo(@FieldType(Vtable, pre ++ name)).optional.child {
                    const Fn = @typeInfo(@FieldType(Vtable, pre ++ name)).optional.child;
                    return @extern(Fn, .{ .name = util.kernelSymbol(@field(Slot, pre ++ name), core), .visibility = .hidden });
                }
                fn unroll(comptime pre: []const u8, comptime _: type, comptime which: []const u8) i32 {
                    return @extern(*const i32, .{ .name = util.kernelSymbol(@field(Slot, pre ++ "gemm_unroll_" ++ which), core), .visibility = .hidden }).*;
                }
                fn unroll3m(comptime pre: []const u8, comptime _: type, comptime which: []const u8) i32 {
                    return @extern(*const i32, .{ .name = util.kernelSymbol(@field(Slot, pre ++ "gemm3m_unroll_" ++ which), core), .visibility = .hidden }).*;
                }
            });
        }
    };
}
