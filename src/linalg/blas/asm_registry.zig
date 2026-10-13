//! The per-core kernels that are assembly upstream, for any routine at any
//! BLAS level. `index` is one architecture's generated `asm/__root__.zig` for
//! a routine: which converted file has which kernel, each file saying which
//! cores (and build variants) it serves. Each kernel is the source's own
//! function, naked, so it is exported and installed as is.

const std = @import("std");
const builtin = @import("builtin");
const util = @import("util.zig");
const Vtable = @import("Vtable.zig");

/// The build variant (`tools/blas_asm/arch/<arch>.zig`) the target needs: the
/// calling convention or pointer width the sources have a path for. It is
/// the compiled target's, so `has`, which the build script asks about the
/// target it is building for, does not depend on it.
pub const variant: []const u8 = switch (builtin.cpu.arch) {
    .x86_64 => if (builtin.os.tag == .windows) "win" else "sysv",
    .arm, .armeb, .thumb, .thumbeb => if (builtin.abi.float() == .hard) "hf" else "sf",
    .powerpc, .powerpcle, .powerpc64, .powerpc64le, .sparc, .sparc64 => if (@sizeOf(usize) == 8) "64" else "32",
    else => "",
};

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

/// Points the vtable at `core`'s kernels from `modules` (an arch's list in
/// `kernels/<arch>.zig`). A core with no linked object uses the build's pick
/// for the target instead, if any; otherwise the portable kernels stay.
pub fn install(vt: *Vtable, comptime modules: anytype, core: anytype) void {
    const Core = @TypeOf(core);
    const pick: Core = if (isLinked(core)) core else if (comptime linked.fallback) |f|
        comptime std.meta.stringToEnum(Core, f).?
    else
        return;
    switch (pick) {
        inline else => |c| {
            if (comptime !isLinked(c)) return;
            inline for (modules) |m| m.install(vt, c);
        },
    }
}

/// Points the vtable slot of each of `roles` that `core` has (`hasRole`) at
/// the kernel the per-core object exports for it, for a routine whose kernels
/// are not all converted assembly (`Kernels.install` does this from `index`).
pub fn installRoles(vt: *Vtable, comptime core: anytype, comptime roles: []const []const u8, comptime hasRole: anytype) void {
    inline for (roles) |role| {
        if (comptime hasRole(core, role)) {
            const Fn = @typeInfo(@FieldType(Vtable, role ++ "_k")).optional.child;
            @field(vt, role ++ "_k") = @extern(Fn, .{ .name = util.kernelSymbol(role, @tagName(core)), .visibility = .hidden });
        }
    }
}

pub fn Kernels(comptime Core: type, comptime index: type) type {
    return struct {
        /// `e`'s function for `core` on this target, if it has one.
        fn kernel(comptime e: anytype, comptime core: Core) ?type {
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
                if (comptime kernel(e, core)) |K| @export(&@field(K, e.role), .{ .name = symbol(e.role, core), .visibility = .hidden });
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
                    const field = e.role ++ "_k";
                    if (comptime @hasField(Vtable, field)) {
                        const Fn = @typeInfo(@FieldType(Vtable, field)).optional.child;
                        @field(vt, field) = @extern(Fn, .{ .name = symbol(e.role, core), .visibility = .hidden });
                    }
                }
            }
        }
    };
}
