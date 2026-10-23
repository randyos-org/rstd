//! The per-core kernels that are assembly upstream, for any routine at any
//! BLAS level. `index` is one architecture's generated `asm/__root__.zig` for
//! a routine: which converted file has which kernel, each file saying which
//! cores (and build variants) it serves. Each kernel is the source's own
//! function, naked, so it is exported and installed as is.

const std = @import("std");
const sysinfo = @import("builtin");
const util = @import("util.zig");
const Vtable = @import("Vtable.zig");
const abi = @import("abi.zig");
const select = @import("cpu/_select.zig");
const vector = @import("vector.zig");

/// The build variant (`tools/blas_asm/arch/<arch>.zig`) the target's calling
/// convention needs (`abi.zig`). The compiled target's, so `has`, which the
/// build script asks about the target it is building for, does not depend on it.
pub const variant: []const u8 = abi.variantName(&sysinfo.target);

/// The cores whose kernel objects this build links, and the one standing in
/// for a detected core that is not linked (`build/blaskern.zig`).
const linked = @import("blaskern_linked");

/// Whether `core`'s kernel object is linked.
pub fn isLinked(core: anytype) bool {
    @setEvalBranchQuota(100_000);
    for (linked.cores) |name| {
        if (std.mem.eql(u8, name, @tagName(core))) return true;
    }
    return false;
}

/// The linked objects as `select.served` reads them: which cores have one, and
/// the build's pick for the target as the last resort.
pub fn LinkedCtx(comptime Core: type) type {
    return struct {
        pub const fallback: ?Core = if (linked.fallback) |f| std.meta.stringToEnum(Core, f) else null;
        pub fn ok(_: @This(), core: Core) bool {
            return isLinked(core);
        }
    };
}

/// The core whose object serves `core`: itself if linked, else the closest linked
/// core that needs nothing `core` does not (`select.served`), else the build's
/// pick for the target. Null means the portable kernels.
pub fn served(core: anytype) ?@TypeOf(core) {
    const Core = @TypeOf(core);
    return select.served(Core, core, LinkedCtx(Core){});
}

/// Points the vtable at `core`'s kernels from `modules` (an arch's list in
/// `kernels/<arch>.zig`), or at those of the core that serves it (`served`).
/// With none, the portable kernels stay.
pub fn install(vt: *Vtable, comptime modules: anytype, core: anytype) void {
    switch (served(core) orelse return) {
        inline else => |c| {
            if (comptime !isLinked(c)) return;
            inline for (modules) |m| m.install(vt, c);
        },
    }
}

/// The vtable slot of `role`: level 1's and `ger`'s are `<role>_k` (`sdot_k`,
/// `cgeru_k`), the other level-2 ones the role itself (`sgemv_n`, `chemv_V`).
fn slot(comptime role: []const u8) []const u8 {
    return if (@hasField(Vtable, role ++ "_k")) role ++ "_k" else role;
}

/// Points the vtable slot of each of `roles` that `core` has (`hasRole`) at
/// the kernel the per-core object exports for it, for a routine whose kernels
/// are not all converted assembly (`Kernels.install` does this from `index`).
pub fn installRoles(vt: *Vtable, comptime core: anytype, comptime roles: []const []const u8, comptime hasRole: anytype) void {
    inline for (roles) |role| {
        if (comptime hasRole(core, role)) {
            const Fn = @typeInfo(@FieldType(Vtable, slot(role))).optional.child;
            @field(vt, slot(role)) = @extern(Fn, .{ .name = util.kernelSymbol(role, @tagName(core)), .visibility = .hidden });
        }
    }
}

pub fn Kernels(comptime Core: type, comptime index: type) type {
    return struct {
        /// `e`'s function for `core` on this target, if it has one and its
        /// role is not `routed` to the portable or `@Vector` kernels.
        fn kernel(comptime e: anytype, comptime core: Core) ?type {
            if (comptime vector.routed(e.role)) return null;
            return e.module.kernel(@tagName(core), variant);
        }

        /// Whether any of these kernels serves `core`, built any way.
        pub fn has(comptime core: Core) bool {
            @setEvalBranchQuota(1_000_000);
            inline for (index.entries) |e| {
                if (comptime e.module.serves(@tagName(core))) return true;
            }
            return false;
        }

        /// Whether `role` has a kernel on `core` for this target.
        pub fn hasRole(comptime core: Core, comptime role: []const u8) bool {
            @setEvalBranchQuota(1_000_000);
            inline for (index.entries) |e| {
                if (comptime std.mem.eql(u8, e.role, role) and kernel(e, core) != null) return true;
            }
            return false;
        }

        /// Whether `role` has a converted kernel on `core` for this target,
        /// `routed` or not.
        pub fn converted(comptime core: Core, comptime role: []const u8) bool {
            @setEvalBranchQuota(1_000_000);
            inline for (index.entries) |e| {
                if (comptime std.mem.eql(u8, e.role, role) and e.module.kernel(@tagName(core), variant) != null) return true;
            }
            return false;
        }

        /// Exported symbol of `role` for `core`.
        pub fn symbol(comptime role: []const u8, comptime core: Core) []const u8 {
            return util.kernelSymbol(role, @tagName(core));
        }

        /// Exports the kernel `role` has on `core`, if it has one.
        pub fn exportOne(comptime role: []const u8, comptime core: Core) void {
            @setEvalBranchQuota(1_000_000);
            inline for (index.entries) |e| {
                if (comptime std.mem.eql(u8, e.role, role)) {
                    if (comptime kernel(e, core)) |K| @export(&@field(K, e.role), .{ .name = symbol(e.role, core), .visibility = .hidden });
                }
            }
        }

        /// Exports every kernel `core` has. Call from a container-level `comptime`
        /// block of the object's root file.
        pub fn exportCore(comptime core: Core) void {
            @setEvalBranchQuota(1_000_000);
            inline for (index.entries) |e| {
                if (comptime kernel(e, core)) |K| {
                    @export(&@field(K, e.role), .{ .name = symbol(e.role, core), .visibility = .hidden });
                } else if (comptime e.module.serves(@tagName(core)) and !hasRole(core, e.role) and !vector.routed(e.role)) {
                    @compileError("no '" ++ variant ++ "' build of " ++ e.role ++ " for " ++ @tagName(core) ++ "; tools/blas_asm/arch and abi.zig disagree");
                }
            }
        }

        /// Points the vtable at the kernels `core` has.
        ///
        /// A role with no matching vtable field (e.g. `sdsdot`: OpenBLAS compiles
        /// it per-core same as `dsdot`, but its own `gotoblas_t` has no
        /// `sdsdot_k` slot either) is still exported, just not installed here.
        pub fn install(vt: *Vtable, comptime core: Core) void {
            @setEvalBranchQuota(1_000_000);
            inline for (index.entries) |e| {
                if (comptime kernel(e, core) != null) {
                    const field = comptime slot(e.role);
                    if (comptime @hasField(Vtable, field)) {
                        const Fn = @typeInfo(@FieldType(Vtable, field)).optional.child;
                        @field(vt, field) = @extern(Fn, .{ .name = symbol(e.role, core), .visibility = .hidden });
                    }
                }
            }
        }
    };
}
