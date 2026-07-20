pub const rstd = @import("../src/__root__.zig");
pub const rstdlib = @import("rstdlib.zig");

const rstdbuild = rstd.buildutils;

pub fn build(b: *rstdbuild.Build) void {
    // standard build setup
    // const optimize = rstdbuild.getOptimize(b);
    const target = rstdbuild.getTarget(b);
    const docs = rstdbuild.addDocs(b);

    // build steps
    const build_options = rstdbuild.addBuildOptionsModule(b);
    _ = rstdlib.addRstd(b, docs, build_options, target);
}
