//! The per-core kernels that are assembly upstream, for any routine at any
//! BLAS level: `index` is that architecture's generated `asm/__root__.zig`
//! (which converted file serves which cores), turned into what the per-core
//! kernel objects (`src/linalg/blas/_kernel.zig`) and the architecture's
//! `install.zig` ask for. Both generics here read `index` the same way; each
//! routine's `<arch>/kernel.zig` picks whichever matches its kernel's calling
//! convention, three directories down from here
//! (`blas1/<routine>/<arch>/kernel.zig`, and the same depth under `blas2/` and
//! `blas3/`).
//!
//! `Kernels` is for a naked function that keeps the source's own ABI, which
//! every routine but `axpy` converts to; its roles come straight off `index`,
//! so it takes nothing else. `Registry` is for a kernel whose operands are
//! pinned to fixed registers by an `asm` block's operand list -- so far only
//! `axpy`'s own kernels, whose signature (and so pin list) is uniform across
//! every core -- and takes the routine's fixed role list (`axpy`'s own copy is
//! `blas1/axpy/kernel.zig`'s `roles`) as a parameter, since nothing here can derive it.

const std = @import("std");
const util = @import("util.zig");
const Vtable = @import("Vtable.zig");

/// `variant` picks between files a source is converted to more than once, such
/// as 32-bit and 64-bit; empty when there is only one.
pub fn Kernels(comptime Core: type, comptime index: type, comptime variant: []const u8) type {
    return struct {
        fn serves(comptime e: anytype, comptime core: Core) bool {
            if (!(e.variant.len == 0 or std.mem.eql(u8, e.variant, variant))) return false;
            inline for (e.cores) |c| {
                if (comptime std.mem.eql(u8, c, @tagName(core))) return true;
            }
            return false;
        }

        /// Whether `core` has any of these kernels.
        pub fn has(comptime core: Core) bool {
            @setEvalBranchQuota(1_000_000);
            inline for (index.entries) |e| {
                if (comptime serves(e, core)) return true;
            }
            return false;
        }

        /// Exports every kernel `core` has. Call from a container-level `comptime`
        /// block of the object's root file.
        pub fn exportCore(comptime core: Core) void {
            @setEvalBranchQuota(1_000_000);
            inline for (index.entries) |e| {
                if (comptime serves(e, core)) {
                    @export(&@field(e.module, e.role), .{ .name = util.kernelSymbol(e.role, @tagName(core)) });
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
                if (comptime serves(e, core)) {
                    const field = e.role ++ "_k";
                    if (comptime @hasField(Vtable, field)) {
                        const Slot = @FieldType(Vtable, field);
                        const Fn = @typeInfo(Slot).optional.child;
                        @field(vt, field) = @extern(Fn, .{ .name = util.kernelSymbol(e.role, @tagName(core)) });
                    }
                }
            }
        }
    };
}

/// `variant` picks between the files a source is converted to more than once,
/// such as 32-bit and 64-bit; empty when there is only one. `roles` is the
/// routine's fixed set of kernel names (`axpy`'s is `saxpy`, `daxpy`, `caxpy`,
/// `caxpyc`, `zaxpy`, `zaxpyc`): unlike `Kernels`, nothing here can read it off
/// `index`, since a pinned-operand kernel's whole-function form carries no
/// entry beyond the file that already names its own role.
pub fn Registry(comptime CoreType: type, comptime index: type, comptime variant: []const u8, comptime roles: []const []const u8) type {
    return struct {
        pub const Core = CoreType;

        /// The file `role` is on `core`, if it is assembly there.
        fn source(comptime core: Core, comptime role: []const u8) ?type {
            return index.source(role, @tagName(core), variant);
        }

        /// Whether `role` has an assembly kernel on `core`.
        pub fn hasRole(comptime core: Core, comptime role: []const u8) bool {
            return comptime source(core, role) != null;
        }

        /// Whether `core` has any assembly kernel of `roles`.
        pub fn hasAny(comptime core: Core) bool {
            inline for (roles) |role| {
                if (comptime source(core, role) != null) return true;
            }
            return false;
        }

        /// Exported symbol of `role` for `core`.
        pub fn symbol(comptime role: []const u8, comptime core: Core) []const u8 {
            return util.kernelSymbol(role, @tagName(core));
        }

        /// Exports the kernel `role` has on `core`, which must have one.
        pub fn exportOne(comptime role: []const u8, comptime core: Core) void {
            @export(&@field(source(core, role).?, role), .{ .name = symbol(role, core) });
        }

        /// Exports every kernel of `roles` that `core` has. Call from a
        /// container-level `comptime` block of the object's root file.
        pub fn exportCore(comptime core: Core) void {
            inline for (roles) |role| {
                if (comptime source(core, role) != null) exportOne(role, core);
            }
        }
    };
}
