//! `Aᵀ y` variants for the `sparse` section (T151): four running sums per
//! column with scalar loads of `x`, and the same with AVX2's `vgatherdpd`
//! (x86_64 with AVX2 only): whether a gather beats scalar loads, apart from
//! the four sums' own gain. Bench-only: the product itself runs four sums
//! now (`sparse/mul.zig`, issue #547), in an order of its own.

const std = @import("std");
const sysinfo = @import("builtin");
const rstd = @import("rstd");

const Csc = rstd.linalg.sparse.Csc(f64);
const V = @Vector(4, f64);

/// Whether this build has the gather variant.
pub const has_gather = sysinfo.cpu.arch == .x86_64 and std.Target.x86.featureSetHas(sysinfo.cpu.features, .avx2);

/// `y_j := A(:, j) · x`, in four sums of every fourth nonzero.
pub fn dots4(a: Csc, x: []const f64, y: []f64) void {
    dots(false, a, x, y);
}

/// `dots4` with `x`'s four values a `vgatherdpd`.
pub fn dotsGather(a: Csc, x: []const f64, y: []f64) void {
    if (comptime !has_gather) unreachable;
    dots(true, a, x, y);
}

fn dots(comptime gathered: bool, a: Csc, x: []const f64, y: []f64) void {
    for (0..a.cols) |j| {
        const k1: usize = a.col_start[j + 1];
        var k: usize = a.col_start[j];
        var acc: V = @splat(0);
        while (k + 4 <= k1) : (k += 4) {
            const r = a.row[k..][0..4];
            const xv: V = if (gathered) gather(x.ptr, r.*) else .{ x[r[0]], x[r[1]], x[r[2]], x[r[3]] };
            acc += @as(V, a.value[k..][0..4].*) * xv;
        }
        var t = @reduce(.Add, acc);
        while (k < k1) : (k += 1) t += a.value[k] * x[a.row[k]];
        y[j] = t;
    }
}

/// `base[idx[i]]` for each lane: a full mask in `ymm15`, consumed by the gather.
inline fn gather(base: [*]const f64, idx: @Vector(4, u32)) V {
    return asm (
        \\vpcmpeqd %%ymm15, %%ymm15, %%ymm15
        \\vxorpd %[r], %[r], %[r]
        \\vgatherdpd %%ymm15, (%[b],%[i],8), %[r]
        : [r] "=&x" (-> V),
        : [b] "r" (base),
          [i] "x" (idx),
        : .{ .ymm15 = true, .memory = true });
}
