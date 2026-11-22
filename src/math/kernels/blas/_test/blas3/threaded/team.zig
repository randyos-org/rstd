//! Threaded level 3 on a caller-owned `Team`'s resident workers, bit for bit.

const zigstd = @import("../../../../../../zigstd.zig");
const threading = @import("../../../threading.zig");
const test_io = @import("../../io.zig");
const types = @import("../../../types.zig");
const th = @import("../threaded.zig");
const trt = @import("routines.zig");

const Depths = th.Depths;
const precisions = th.precisions;
const shapes = th.shapes;
const routineShapes = trt.routineShapes;
const routines = trt.routines;
const shrunkCore = trt.shrunkCore;

test "a Team's level 3 keeps the serial bits: gemm and the other routines, on resident workers" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    var threaded: test_io.TestIo = .init(zigstd.testAllocator, .{ .async_limit = .limited(8) });
    defer threaded.deinit();
    var team: threading.Team = undefined;
    team.init(threaded.io(), 4);
    defer team.deinit();
    try zigstd.testExpectTrue(team.size() > 1);
    inline for (precisions) |T| {
        try shapes(T, team.io());
        const vt = Depths.spied(T, shrunkCore(T));
        for (routines) |r| {
            if (!types.isComplex(T) and (r == .hemm or r == .herk or r == .her2k)) continue;
            try routineShapes(T, &vt, team.io(), r);
        }
    }
}
