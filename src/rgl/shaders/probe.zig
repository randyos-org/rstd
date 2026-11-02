//! One fragment shader per shaderlib math operation, so each one's SPIR-V
//! can be disassembled on its own (`zig build shader-probes -Dvulkan`;
//! lib/readme.md). The root of every probe; `probe_op` names the `Op`.
//! The build imports this file for `Op` alone.

const std = @import("std");
const sysinfo = @import("builtin");
const spirv = std.spirv;
const rstd = @import("rstd");
const shader_mod = rstd.rgl.shaders;
const lib = shader_mod.lib;
const types = shader_mod.types;
const math = lib.math;
const Vec3 = types.Vec3;
const Vec4 = types.Vec4;
const Mat4 = types.Mat4;
const param = lib.compile.param;

pub const shader_opts = @import("shader_opts");

comptime {
    // not when the build script imports this for `Op`
    if (sysinfo.cpu.arch.isSpirV()) _ = Probe(@field(Op, @import("probe_op").op)).main;
}

pub const Op = enum {
    dot,
    dot_reduce,
    dot_asm,
    cross,
    cross_asm,
    length,
    distance,
    normalize,
    reflect,
    refract,
    faceforward,
    inversesqrt,
    inversesqrt_vec,
    fma_muladd,
    fma_builtin,
    mix,
    mix_scalar,
    step,
    smoothstep,
    fract,
    fract_scalar,
    clamp,
    clamp_vec,
    minmax,
    mat_vec,
    vec_mat,
    mat_mat,
    vec_arith,
};

const in_a = param("a", Vec3, 0);
const in_b = param("b", Vec3, 1);
const in_c = param("c", Vec3, 2);
const in_t = param("t", f32, 3);
const in_m0 = param("m0", Vec4, 4);
const in_m1 = param("m1", Vec4, 5);
const in_m2 = param("m2", Vec4, 6);
const in_m3 = param("m3", Vec4, 7);
const in_v = param("v", Vec4, 8);

fn out3(v: Vec3) Vec4 {
    return .{ v[0], v[1], v[2], 1 };
}

fn out1(s: f32) Vec4 {
    return @splat(s);
}

/// GLSL.std.450 Cross (68) through inline asm. A literal template:
/// `std.fmt.comptimePrint` here crashes the SPIR-V backend (S1).
fn crossAsm(a: Vec3, b: Vec3) Vec3 {
    return asm (
        \\%set = OpExtInstImport "GLSL.std.450"
        \\%ret = OpExtInst %T %set 68 %a %b
        : [ret] "" (-> Vec3),
        : [T] "t" (Vec3),
          [a] "" (a),
          [b] "" (b),
    );
}

fn dotAsm(a: Vec3, b: Vec3) f32 {
    return asm (
        \\%ret = OpDot %f32 %a %b
        : [ret] "" (-> f32),
        : [f32] "t" (f32),
          [a] "" (a),
          [b] "" (b),
    );
}

/// GLSL's `m * v` on a column-major matrix.
fn matVec(m: Mat4, v: Vec4) Vec4 {
    var r = m[0] * @as(Vec4, @splat(v[0]));
    inline for (1..4) |i| r += m[i] * @as(Vec4, @splat(v[i]));
    return r;
}

/// GLSL's `v * m`: one dot per column.
fn vecMat(v: Vec4, m: Mat4) Vec4 {
    var r: Vec4 = undefined;
    inline for (0..4) |i| r[i] = math.dot(v, m[i]);
    return r;
}

fn run(comptime op: Op) Vec4 {
    const a = in_a.*;
    const b = in_b.*;
    const c = in_c.*;
    const t = in_t.*;
    const m: Mat4 = .{ in_m0.*, in_m1.*, in_m2.*, in_m3.* };
    const v = in_v.*;
    return switch (op) {
        .dot => out1(math.dot(a, b)),
        .dot_reduce => out1(@reduce(.Add, a * b)),
        .dot_asm => out1(dotAsm(a, b)),
        .cross => out3(math.cross(a, b)),
        .cross_asm => out3(crossAsm(a, b)),
        .length => out1(math.length(a)),
        .distance => out1(math.distance(a, b)),
        .normalize => out3(math.normalize(a)),
        .reflect => out3(math.reflect(a, b)),
        .refract => out3(math.refract(a, b, t)),
        .faceforward => out3(math.faceforward(a, b, c)),
        .inversesqrt => out1(math.inversesqrt(t)),
        .inversesqrt_vec => out3(math.inversesqrt(a)),
        .fma_muladd => out3(rstd.linalg.vector.mulAdd(Vec3, a, b, c)),
        .fma_builtin => out3(@mulAdd(Vec3, a, b, c)),
        .mix => out3(math.mix(a, b, t)),
        .mix_scalar => out1(math.mix(a[0], b[0], t)),
        .step => out1(math.step_scalar(a[0], t)),
        .smoothstep => out1(math.smoothstep_scalar(a[0], b[0], t)),
        .fract => out3(math.fract(a)),
        .fract_scalar => out1(math.fract(t)),
        .clamp => out1(math.clamp(t, 0.0, 1.0)),
        .clamp_vec => out3(math.clamp(a, b, c)),
        .minmax => out3(@max(@min(a, b), c)),
        .mat_vec => matVec(m, v),
        .vec_mat => vecMat(v, m),
        .mat_mat => blk: {
            const n: Mat4 = .{ v, in_m1.*, in_m2.*, in_m0.* };
            var r: Mat4 = undefined;
            inline for (0..4) |j| r[j] = matVec(m, n[j]);
            break :blk r[0] + r[1] + r[2] + r[3];
        },
        .vec_arith => out3(a * b + c),
    };
}

pub fn Probe(comptime op: Op) type {
    return struct {
        fn vertMain() callconv(.spirv_vertex) void {
            spirv.position_out.* = @splat(0);
        }

        fn fragMain() callconv(.{ .spirv_fragment = .{} }) void {
            lib.compile.frag_color_out.* = run(op);
        }

        pub const main = lib.compile.exportShaderMain(struct {
            pub const vert = vertMain;
            pub const frag = fragMain;
        });
    };
}
