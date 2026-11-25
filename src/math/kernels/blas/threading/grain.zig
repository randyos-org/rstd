//! The work a team needs per worker (`threading.md`).

const sysinfo = @import("builtin");
const types = @import("../types.zig");

/// Real multiply-adds (a complex one counts four) a level-3 call needs per
/// worker, public (#591) or inside a threaded routine (#546). Measured, not
/// upstream's (`threading.md`): a team costs 2 us a worker to start on the
/// M6 (macOS), where `potrf`'s rate is flat from 400000 to 2500000, and 4 to
/// 10 on Zen 3 (Windows), where teams of two lose at 1000000; upstream's
/// workers spin between calls. Elsewhere macOS's, unmeasured. No result
/// depends on it.
pub const team_grain: u64 = if (sysinfo.target.os.tag == .windows) 4_000_000 else 1_000_000;

/// Bytes of matrix a level-2 call (`ger`, `symv`, `hemv`) needs per worker.
/// Measured on the M6, not upstream's: memory bound, they peak at 4 to 8
/// workers at 4 MiB and lose past them, a worker costing its start (#707).
pub const matrix_grain: u64 = 640 * 1024;

/// The workers a level-3 call of `madds` real multiply-adds (a complex one
/// counting four) takes: `want`, at most one per `grain`. The drivers give
/// the same bits for any count, so this changes none.
pub fn crew(comptime T: type, want: u32, grain: u64, madds: u64) u32 {
    const each: u64 = if (comptime types.isComplex(T)) 4 else 1;
    return share(want, grain, madds *| each);
}

/// `want` workers, at most one per `grain` of `work`, at least one.
pub fn share(want: u32, grain: u64, work: u64) u32 {
    return @intCast(@min(want, @max(1, work / @max(grain, 1))));
}
