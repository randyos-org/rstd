// const std = @import("std");
const rstdbuild = @import("__root__.zig");
const Build = rstdbuild.Build;
const Module = rstdbuild.Module;
const ResolvedTarget = rstdbuild.ResolvedTarget;

/// Shallow-clones `existing` into a new `Module` bound to `target`.
///
/// A `Compile` step needs its root module to have a known target, but a
/// module meant for `addImport` (like `dep.module(...)`) is usually
/// target-less so each consumer picks its own. This reuses `existing`'s
/// contents for a standalone `Compile` step (docs, tests) without
/// redeclaring root_source_file/imports/options.
pub fn cloneModuleWithTarget(b: *Build, existing: *Module, target: ?ResolvedTarget) !*Module {
    const clone = try b.allocator.create(Module);
    clone.* = existing.*;
    clone.resolved_target = target;
    // The struct copy above also copies the memoized import-graph cache
    // (keyed off existing's identity as root); reset it so the clone -
    // a distinct root - recomputes its own instead of reusing existing's.
    clone.cached_graph = .{ .modules = &.{}, .names = &.{} };
    // Also deep-copy `import_table`: a plain struct copy shares its
    // backing array with `existing`, so a later `addImport` on the clone
    // and one on `existing` (or another clone of it) can land in the same
    // slot and silently overwrite each other once either grows in place.
    clone.import_table = try existing.import_table.clone(b.allocator);
    return clone;
}
