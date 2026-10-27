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

const Slot = util.Slot;

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

/// One row of a generated `asm/__root__.zig`: a converted file, and the slot
/// its kernel fills.
pub const Entry = struct { slot: Slot, module: type };

/// Exports `f` as `core`'s kernel for `slot`, under `util.kernelSymbol`.
pub fn exportKernel(comptime f: anytype, comptime slot: Slot, comptime core: anytype) void {
    @export(f, .{ .name = util.kernelSymbol(slot, core), .visibility = .hidden });
}

/// The kernel `core`'s object exports for `slot` (`exportKernel`).
pub fn externKernel(comptime slot: Slot, comptime core: anytype) util.SlotFn(slot) {
    return @extern(util.SlotFn(slot), .{ .name = util.kernelSymbol(slot, core), .visibility = .hidden });
}

/// Points each of `slots` that `core` has (`hasSlot`) at the kernel the
/// per-core object exports for it, for a routine whose kernels are not all
/// converted assembly (`Kernels.install` does this from `index`).
pub fn installSlots(vt: *Vtable, comptime core: anytype, comptime slots: []const Slot, comptime hasSlot: anytype) void {
    inline for (slots) |slot| {
        if (comptime hasSlot(core, slot)) @field(vt, @tagName(slot)) = externKernel(slot, core);
    }
}

pub fn Kernels(comptime Core: type, comptime index: type) type {
    return struct {
        /// `e`'s function for `core` on this target, if it has one and its
        /// slot is not `routed` to the portable or `@Vector` kernels.
        fn kernel(comptime e: Entry, comptime core: Core) ?type {
            if (comptime vector.routed(e.slot)) return null;
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

        /// Whether `slot` has a kernel on `core` for this target.
        pub fn hasSlot(comptime core: Core, comptime slot: Slot) bool {
            inline for (index.entries) |e| {
                if (comptime e.slot == slot and kernel(e, core) != null) return true;
            }
            return false;
        }

        /// Whether `slot` has a converted kernel on `core` for this target,
        /// `routed` or not.
        pub fn converted(comptime core: Core, comptime slot: Slot) bool {
            inline for (index.entries) |e| {
                if (comptime e.slot == slot and e.module.kernel(@tagName(core), variant) != null) return true;
            }
            return false;
        }

        /// Exports the kernel `slot` has on `core`, if it has one.
        pub fn exportOne(comptime slot: Slot, comptime core: Core) void {
            inline for (index.entries) |e| {
                if (comptime e.slot == slot) {
                    if (comptime kernel(e, core)) |K| exportKernel(&K.entry, e.slot, core);
                }
            }
        }

        /// Exports every kernel `core` has. Call from a container-level `comptime`
        /// block of the object's root file.
        pub fn exportCore(comptime core: Core) void {
            inline for (index.entries) |e| {
                if (comptime kernel(e, core)) |K| {
                    exportKernel(&K.entry, e.slot, core);
                } else if (comptime e.module.serves(@tagName(core)) and !hasSlot(core, e.slot) and !vector.routed(e.slot)) {
                    @compileError("no '" ++ variant ++ "' build of " ++ @tagName(e.slot) ++ " for " ++ @tagName(core) ++ "; tools/blas_asm/arch and abi.zig disagree");
                }
            }
        }

        /// Points the vtable at the kernels `core` has.
        pub fn install(vt: *Vtable, comptime core: Core) void {
            inline for (index.entries) |e| {
                if (comptime kernel(e, core) != null) @field(vt, @tagName(e.slot)) = externKernel(e.slot, core);
            }
        }
    };
}
