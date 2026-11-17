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
const matrix = rstd.linalg.matrix;
const Vec3 = types.Vec3;
const Vec4 = types.Vec4;
const Mat4 = types.Mat4;
const param = lib.compile.paramFlat;
const input = lib.compile.input;

pub const shader_opts = @import("shader_opts");

comptime {
    // not when the build script imports this for `Op`
    if (sysinfo.cpu.arch.isSpirv()) _ = Probe(@field(Op, @import("probe_op").op)).main;
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
    det,
    inv,
    outer,
    dot2,
    vec_arith,
};

/// What `spirv-dis` must show in `op`'s module: its native instruction, an
/// extended one between spaces, or null where none is due (lib/readme.md's
/// table). `shader-probes` fails on a probe that lacks it.
pub fn expect(op: Op) ?[]const u8 {
    return switch (op) {
        .dot, .dot_asm, .dot2 => "OpDot",
        .cross, .cross_asm => " Cross ",
        .length => " Length ",
        .distance => " Distance ",
        .normalize => " Normalize ",
        .reflect => " Reflect ",
        .refract => " Refract ",
        .faceforward => " FaceForward ",
        .inversesqrt, .inversesqrt_vec => " InverseSqrt ",
        .fma_muladd, .fma_builtin => " Fma ",
        .mix, .mix_scalar => " FMix ",
        .step => " Step ",
        .smoothstep => " SmoothStep ",
        .fract, .fract_scalar => " Fract ",
        .clamp, .clamp_vec => " FClamp ",
        .minmax => " FMin ",
        .mat_vec => "OpMatrixTimesVector",
        .vec_mat => "OpVectorTimesMatrix",
        .mat_mat => "OpMatrixTimesMatrix",
        .det => " Determinant ",
        .inv => " MatrixInverse ",
        .outer => "OpOuterProduct",
        .dot_reduce, .vec_arith => null,
    };
}

// flat, so the fragment stage reads the instance's values exactly; the
// harness (rgl/vulkan/_test) feeds them as instance attributes 0 to 8
const in_a = param("a", Vec3, 0);
const in_b = param("b", Vec3, 1);
const in_c = param("c", Vec3, 2);
const in_t = param("t", f32, 3);
const in_m0 = param("m0", Vec4, 4);
const in_m1 = param("m1", Vec4, 5);
const in_m2 = param("m2", Vec4, 6);
const in_m3 = param("m3", Vec4, 7);
const in_v = param("v", Vec4, 8);

const attr_a = input("attr_a", Vec3, 0);
const attr_b = input("attr_b", Vec3, 1);
const attr_c = input("attr_c", Vec3, 2);
const attr_t = input("attr_t", f32, 3);
const attr_m0 = input("attr_m0", Vec4, 4);
const attr_m1 = input("attr_m1", Vec4, 5);
const attr_m2 = input("attr_m2", Vec4, 6);
const attr_m3 = input("attr_m3", Vec4, 7);
const attr_v = input("attr_v", Vec4, 8);

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

/// A probe's inputs: the vertex attributes, as the harness feeds them.
pub const Inputs = struct {
    a: Vec3,
    b: Vec3,
    c: Vec3,
    t: f32,
    m: [4]Vec4,
    v: Vec4,
};

fn run(comptime op: Op) Vec4 {
    return compute(op, .{
        .a = in_a.*,
        .b = in_b.*,
        .c = in_c.*,
        .t = in_t.*,
        .m = .{ in_m0.*, in_m1.*, in_m2.*, in_m3.* },
        .v = in_v.*,
    });
}

/// `op` on `in`: the shader's result, or on the host the CPU's, against which
/// the harness checks the shader's (`_asm` ops by their portable forms).
pub fn compute(comptime op: Op, in: Inputs) Vec4 {
    const on_spirv = comptime sysinfo.cpu.arch.isSpirv();
    const a = in.a;
    const b = in.b;
    const c = in.c;
    const t = in.t;
    const m = matrix.fromCols(Mat4, .{ in.m[0], in.m[1], in.m[2], in.m[3] });
    const v = in.v;
    return switch (op) {
        .dot => out1(math.dot(a, b)),
        .dot_reduce => out1(@reduce(.Add, a * b)),
        .dot_asm => out1(if (on_spirv) dotAsm(a, b) else math.dot(a, b)),
        .cross => out3(math.cross(a, b)),
        .cross_asm => out3(if (on_spirv) crossAsm(a, b) else math.cross(a, b)),
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
        .step => out1(math.stepScalar(a[0], t)),
        .smoothstep => out1(math.smoothstepScalar(a[0], b[0], t)),
        .fract => out3(math.fract(a)),
        .fract_scalar => out1(math.fract(t)),
        .clamp => out1(math.clamp(t, 0.0, 1.0)),
        .clamp_vec => out3(math.clamp(a, b, c)),
        .minmax => out3(@max(@min(a, b), c)),
        .mat_vec => matrix.mul(m, v),
        .vec_mat => matrix.mul(v, m),
        .mat_mat => blk: {
            const r = matrix.mul(m, matrix.fromCols(Mat4, .{ v, in.m[1], in.m[2], in.m[0] })).data;
            break :blk @as(Vec4, r[0]) + r[1] + r[2] + r[3];
        },
        .det => out1(matrix.det(m)),
        .inv => blk: {
            var r = matrix.identity(Mat4).data;
            if (matrix.inv(m)) |x| r = x.data else |_| {}
            break :blk @as(Vec4, r[0]) + r[1] + r[2] + r[3];
        },
        .outer => blk: {
            const r = matrix.outer(v, in.m[0]).data;
            break :blk @as(Vec4, r[0]) + r[1] + r[2] + r[3];
        },
        .dot2 => out1(math.dot2(a)),
        .vec_arith => out3(a * b + c),
    };
}

pub fn Probe(comptime op: Op) type {
    return struct {
        /// A triangle over the whole target, its inputs passed through.
        fn vertMain() callconv(.spirv_vertex) void {
            const i = spirv.vertex_index;
            const x: f32 = if (i == 1) 3 else -1;
            const y: f32 = if (i == 2) 3 else -1;
            spirv.position_out.* = .{ x, y, 0, 1 };
            in_a.* = attr_a.*;
            in_b.* = attr_b.*;
            in_c.* = attr_c.*;
            in_t.* = attr_t.*;
            in_m0.* = attr_m0.*;
            in_m1.* = attr_m1.*;
            in_m2.* = attr_m2.*;
            in_m3.* = attr_m3.*;
            in_v.* = attr_v.*;
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
