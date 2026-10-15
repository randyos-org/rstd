//! The level-2 interface: the checks, `beta` scaling and pointer set-up
//! OpenBLAS's `interface/{gemv,zgemv,ger,zger,symv,zsymv,zhemv}.c` do before
//! calling a kernel. Like `blas1/interface.zig`, each takes the element type
//! `T` explicitly, with `m`, `n`, `lda` and the strides `BlasInt`; matrices
//! are column-major.
//!
//! A negative stride names the vector's last element first: the interface
//! moves the pointer to it. The kernels need a scratch `buffer`; the
//! interface takes it from a small stack array (OpenBLAS's `STACK_ALLOC`) and
//! only allocates beyond that, which is why these take an allocator.
//!
//! Invalid arguments return an `ArgError` naming the argument, where the
//! reference calls `xerbla` with its number (checked in the same order, so
//! the first bad argument wins). Single-threaded: the reference's threaded
//! drivers would be chosen where marked.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const Vtable = @import("../Vtable.zig");
const check = @import("../blas1/check.zig");
const sym = @import("symmetric.zig");
const portable = @import("../portable.zig");
const guard = @import("guard.zig");

const Allocator = std.mem.Allocator;
const BlasInt = types.BlasInt;
const BlasLong = types.BlasLong;
const XDouble = types.XDouble;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const ComplexXD = types.ComplexXD;

/// The reference's `xerbla` argument, by name.
/// Bad arguments, shared by every level-2 interface file; there is no xerbla.
pub const ArgError = error{
    InvalidUplo,
    InvalidTrans,
    InvalidDiag,
    InvalidM,
    InvalidN,
    InvalidK,
    InvalidLda,
    InvalidIncx,
    InvalidIncy,
};
pub const Error = ArgError || Allocator.Error;

/// `gemv`'s `TRANS`: `n` is `A`, `t` is `A^T`, `r` is `conj(A)`, `c` is
/// `A^H`; complex also has `o` (`A conj(x)`), `u` (`A^T conj(x)`), `s`
/// (`conj(A) conj(x)`) and `d` (`A^H conj(x)`). Real `r` and `c` are `n` and
/// `t`.
pub const Trans = enum { n, t, r, c, o, u, s, d };

/// `symv` and `hemv`'s `UPLO`: the stored triangle. `v` and `m` are
/// `hemv`'s upper and lower triangles of the transposed matrix (`HEMVREV`).
pub const Uplo = enum { u, l, v, m };

fn transposed(t: Trans) bool {
    return switch (t) {
        .t, .c, .u, .d => true,
        .n, .r, .o, .s => false,
    };
}

fn isComplex(comptime T: type) bool {
    return check.isComplex(T);
}

const Float = check.FloatOf;

/// The reference's precision letter for `T` (`f64` before `XDouble`, which
/// may be the same type).
fn letter(comptime T: type) []const u8 {
    if (T == f32 or T == Complex32) return if (T == f32) "s" else "c";
    if (T == f64 or T == Complex64) return if (T == f64) "d" else "z";
    if (T == XDouble or T == ComplexXD) return if (T == XDouble) "q" else "x";
    @compileError("not a BLAS element type: " ++ @typeName(T));
}

/// The kernel in the slot named `letter(T) ++ name`.
fn slot(comptime T: type, comptime name: []const u8) @typeInfo(@FieldType(Vtable, letter(T) ++ name)).optional.child {
    return blas.kernel(@field(blas.Slot, letter(T) ++ name));
}

fn raw(comptime T: type, p: [*]T) [*]Float(T) {
    return @ptrCast(p);
}

fn long(v: BlasInt) BlasLong {
    return @intCast(v);
}

/// Where a vector of `n` with stride `inc` starts: its last element when
/// `inc < 0`.
fn first(comptime T: type, p: [*]T, n: BlasInt, inc: BlasInt) [*]T {
    return if (inc < 0) @ptrCast(util.at(T, p, -(long(n) - 1) * inc)) else p;
}

fn isOne(comptime T: type, v: T) bool {
    return if (comptime isComplex(T)) v.re == 1 and v.im == 0 else v == 1;
}

fn isZero(comptime T: type, v: T) bool {
    return if (comptime isComplex(T)) v.re == 0 and v.im == 0 else v == 0;
}

/// `y := beta*y` over `n` elements at the unadjusted `y`, the way the
/// reference does it: with the scal kernel's flag 0, so a zero `beta`
/// overwrites NaN.
fn scaleY(comptime T: type, n: BlasInt, beta: T, y: [*]T, incy: BlasInt) void {
    const py = raw(T, y);
    const inc: BlasLong = @intCast(@abs(incy));
    if (comptime isComplex(T)) {
        _ = slot(T, "scal_k")(long(n), 0, 0, beta.re, beta.im, py, inc, py, 0, py, 0);
    } else {
        _ = slot(T, "scal_k")(long(n), 0, 0, beta, py, inc, py, 0, py, 0);
    }
}

// ---- scratch buffer ----------------------------------------------------------

/// `MAX_STACK_ALLOC`: up to this many bytes of scratch live on the stack, in
/// every level-2 driver; the rest comes from the allocator.
pub const stack_bytes = 2048;
/// The reference's stack buffers are `aligned(0x20)`; no ported kernel needs
/// more than 16 (`buffers.md`).
pub const stack_align = 32;
/// `blas_memory_alloc`'s buffers are pages; the generic symv and hemv lay
/// their copies out on page boundaries from it.
pub const heap_align = 4096;

/// Scratch for `F` elements: on the stack when `stack_bytes` holds them, else
/// from the allocator.
fn Scratch(comptime F: type) type {
    return struct {
        const Self = @This();
        ptr: [*]F,
        len: usize,
        heap: ?[]align(heap_align) F = null,

        fn init(allocator: Allocator, stack: *align(stack_align) [stack_bytes]u8, len: usize) Allocator.Error!Self {
            if (len * @sizeOf(F) <= stack_bytes) return .{ .ptr = @ptrCast(stack), .len = len };
            const h = try allocator.alignedAlloc(F, .fromByteUnits(heap_align), len);
            return .{ .ptr = h.ptr, .len = len, .heap = h };
        }

        /// Always from the allocator, as the reference's `blas_memory_alloc`.
        fn heapOnly(allocator: Allocator, len: usize) Allocator.Error!Self {
            const h = try allocator.alignedAlloc(F, .fromByteUnits(heap_align), len);
            return .{ .ptr = h.ptr, .len = len, .heap = h };
        }

        fn slice(self: Self) []F {
            return self.ptr[0..self.len];
        }

        fn deinit(self: Self, allocator: Allocator) void {
            if (self.heap) |h| allocator.free(h);
        }
    };
}

// ---- gemv --------------------------------------------------------------------

/// `?GEMV` (`interface/gemv.c`, `zgemv.c`): `y := alpha*op(A)*x + beta*y`,
/// `A` being `m` by `n`. `y` is scaled by `beta` first (a zero `beta`
/// overwrites NaN), and nothing more happens for a zero `alpha`.
pub fn gemv(allocator: Allocator, comptime T: type, trans: Trans, m: BlasInt, n: BlasInt, alpha: T, a: [*]T, lda: BlasInt, x: [*]T, incx: BlasInt, beta: T, y: [*]T, incy: BlasInt) Error!void {
    if (!isComplex(T) and @backingInt(trans) > @backingInt(Trans.c)) return error.InvalidTrans;
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    if (m == 0 or n == 0) return;
    const tr = transposed(trans);
    const lenx = if (tr) m else n;
    const leny = if (tr) n else m;
    if (!isOne(T, beta)) scaleY(T, leny, beta, y, incy);
    if (isZero(T, alpha)) return;
    switch (trans) {
        inline else => |t| {
            // `r`..`d` have no real kernel; `gemv` returned above for those
            if (comptime !isComplex(T) and @backingInt(t) > @backingInt(Trans.c)) unreachable;
            try gemvRun(allocator, T, t, m, n, alpha, a, lda, first(T, x, lenx, incx), incx, first(T, y, leny, incy), incy);
        },
    }
}

fn gemvRun(allocator: Allocator, comptime T: type, comptime t: Trans, m: BlasInt, n: BlasInt, alpha: T, a: [*]T, lda: BlasInt, x: [*]T, incx: BlasInt, y: [*]T, incy: BlasInt) Allocator.Error!void {
    const F = Float(T);
    const complex = comptime isComplex(T);
    if (comptime !complex and @backingInt(t) > @backingInt(Trans.c)) return;
    const name = comptime "gemv_" ++ if (complex) @tagName(t) else if (transposed(t)) "t" else "n";
    var stack: [stack_bytes]u8 align(stack_align) = undefined;
    const sc = try Scratch(F).init(allocator, &stack, sym.gemvLen(F, if (complex) 2 else 1, @intCast(m), @intCast(n)));
    defer sc.deinit(allocator);
    // `cgemv_t.S` returns NaN if the buffer holds NaN or infinity; the
    // reference clears its stack buffer for every `trans` but `n`, here the
    // heap one too since the allocator hands out whatever it holds
    if (comptime complex and t != .n and sysinfo.cpu.arch == .x86_64) @memset(sc.slice(), 0);
    // the reference threads big problems here (`gemv_thread_*`)
    if (comptime complex) {
        _ = slot(T, name)(long(m), long(n), 0, alpha.re, alpha.im, raw(T, a), long(lda), raw(T, x), incx, raw(T, y), incy, sc.ptr);
    } else {
        _ = slot(T, name)(long(m), long(n), 0, alpha, a, long(lda), x, incx, y, incy, sc.ptr);
    }
}

// ---- half-precision gemv -----------------------------------------------------

/// The half-precision `gemv`s: `sbgemv` is bfloat16 in and `f32` out, `bgemv`
/// bfloat16 throughout, `shgemv` `f16` in and `f32` out.
pub const HalfGemv = enum {
    sbgemv,
    bgemv,
    shgemv,

    fn In(comptime k: HalfGemv) type {
        return if (k == .shgemv) types.HF16 else types.BF16;
    }

    /// The type of `y`, `alpha` and `beta`.
    fn Out(comptime k: HalfGemv) type {
        return if (k == .bgemv) types.BF16 else f32;
    }
};

/// `SBGEMV`, `BGEMV`, `SHGEMV` (`interface/sbgemv.c`): `y := alpha*op(A)*x +
/// beta*y`, `A` being `m` by `n`; `trans` is `n`/`r` or `t`/`c`. Unlike
/// `gemv`, the kernel applies `beta` itself, and a zero `alpha` only scales
/// `y` (when `beta` is not 1) without reading `A` or `x`. The kernels take no
/// scratch buffer, so there is no allocator; the reference's threaded
/// drivers are not here.
pub fn gemvHalf(comptime kind: HalfGemv, trans: Trans, m: BlasInt, n: BlasInt, alpha: kind.Out(), a: [*]kind.In(), lda: BlasInt, x: [*]kind.In(), incx: BlasInt, beta: kind.Out(), y: [*]kind.Out(), incy: BlasInt) ArgError!void {
    const In = kind.In();
    const Out = kind.Out();
    if (@backingInt(trans) > @backingInt(Trans.c)) return error.InvalidTrans;
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    if (m == 0 or n == 0) return;
    const tr = transposed(trans);
    const lenx = if (tr) m else n;
    const leny = if (tr) n else m;
    // `SBF16TOS_K` of the scalars, for bfloat16
    const alpha_f = halfToF32(Out, alpha);
    const beta_f = halfToF32(Out, beta);
    if (alpha_f == 0) {
        if (beta_f != 1) {
            const inc: BlasLong = @intCast(@abs(incy));
            if (comptime kind == .bgemv) {
                _ = blas.kernel(.bscal_k)(long(leny), 0, 0, beta, y, inc, y, 0, y, 0);
            } else {
                _ = blas.kernel(.sscal_k)(long(leny), 0, 0, beta, y, inc, y, 0, y, 0);
            }
        }
        return;
    }
    const px = first(In, x, lenx, incx);
    const py = first(Out, y, leny, incy);
    const name = @tagName(kind);
    _ = if (tr)
        blas.kernel(@field(blas.Slot, name ++ "_t"))(long(m), long(n), alpha, a, long(lda), px, incx, beta, py, incy)
    else
        blas.kernel(@field(blas.Slot, name ++ "_n"))(long(m), long(n), alpha, a, long(lda), px, incx, beta, py, incy);
}

/// A `gemvHalf` scalar as `f32`.
fn halfToF32(comptime T: type, v: T) f32 {
    return if (T == f32) v else @import("../ext/bf16/__root__.zig").scalarToF32(v);
}

// ---- ger ---------------------------------------------------------------------

/// `?GER`, `?GERU` (`interface/ger.c`, `zger.c`): `A := alpha*x*y^T + A`.
pub fn ger(allocator: Allocator, comptime T: type, m: BlasInt, n: BlasInt, alpha: T, x: [*]T, incx: BlasInt, y: [*]T, incy: BlasInt, a: [*]T, lda: BlasInt) Error!void {
    return gerImpl(allocator, T, false, m, n, alpha, x, incx, y, incy, a, lda);
}

/// `?GERC`: `A := alpha*x*y^H + A`. Complex only.
pub fn gerc(allocator: Allocator, comptime T: type, m: BlasInt, n: BlasInt, alpha: T, x: [*]T, incx: BlasInt, y: [*]T, incy: BlasInt, a: [*]T, lda: BlasInt) Error!void {
    if (comptime !isComplex(T)) @compileError("gerc: " ++ @typeName(T) ++ " is not complex");
    return gerImpl(allocator, T, true, m, n, alpha, x, incx, y, incy, a, lda);
}

fn gerImpl(allocator: Allocator, comptime T: type, comptime conj: bool, m: BlasInt, n: BlasInt, alpha: T, x: [*]T, incx: BlasInt, y: [*]T, incy: BlasInt, a: [*]T, lda: BlasInt) Error!void {
    const F = Float(T);
    const complex = comptime isComplex(T);
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (m == 0 or n == 0) return;
    if (isZero(T, alpha)) return;
    const name = comptime if (complex) (if (conj) "gerc_k" else "geru_k") else "ger_k";
    const px = first(T, x, m, incx);
    const py = first(T, y, n, incy);
    var stack: [stack_bytes]u8 align(stack_align) = undefined;
    // the kernels use the buffer only to copy a strided `x` (`buffers.md`);
    // for a unit one the reference may pass NULL, so no kernel reads it
    const unit = incx == 1;
    const sc: Scratch(F) = if (unit) .{ .ptr = @ptrCast(&stack), .len = 0 } else try Scratch(F).init(allocator, &stack, (if (complex) 2 else 1) * @as(usize, @intCast(m)));
    defer sc.deinit(allocator);
    // the reference threads big problems here (`ger_thread`)
    if (comptime complex) {
        _ = slot(T, name)(long(m), long(n), 0, alpha.re, alpha.im, raw(T, px), incx, raw(T, py), incy, raw(T, a), long(lda), sc.ptr);
    } else {
        _ = slot(T, name)(long(m), long(n), 0, alpha, px, incx, py, incy, a, long(lda), sc.ptr);
    }
}

// ---- symv, hemv --------------------------------------------------------------

/// `?SYMV` (`interface/symv.c`, `zsymv.c`): `y := alpha*S*x + beta*y` for the
/// symmetric `S` whose `uplo` triangle (`u` or `l`) is stored in `a`.
pub fn symv(allocator: Allocator, comptime T: type, uplo: Uplo, n: BlasInt, alpha: T, a: [*]T, lda: BlasInt, x: [*]T, incx: BlasInt, beta: T, y: [*]T, incy: BlasInt) Error!void {
    if (uplo == .v or uplo == .m) return error.InvalidUplo;
    return mv(allocator, T, "symv_", uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}

/// `?HEMV` (`interface/zhemv.c`): the same for the Hermitian `S`; `v` and
/// `m` use the upper and lower triangles for `conj(S)`. Complex only.
pub fn hemv(allocator: Allocator, comptime T: type, uplo: Uplo, n: BlasInt, alpha: T, a: [*]T, lda: BlasInt, x: [*]T, incx: BlasInt, beta: T, y: [*]T, incy: BlasInt) Error!void {
    if (comptime !isComplex(T)) @compileError("hemv: " ++ @typeName(T) ++ " is not complex");
    return mv(allocator, T, "hemv_", uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}

fn mv(allocator: Allocator, comptime T: type, comptime base: []const u8, uplo: Uplo, n: BlasInt, alpha: T, a: [*]T, lda: BlasInt, x: [*]T, incx: BlasInt, beta: T, y: [*]T, incy: BlasInt) Error!void {
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    if (n == 0) return;
    if (!isOne(T, beta)) scaleY(T, n, beta, y, incy);
    if (isZero(T, alpha)) return;
    const px = first(T, x, n, incx);
    const py = first(T, y, n, incy);
    switch (uplo) {
        inline else => |u| try mvRun(allocator, T, base ++ switch (u) {
            .u => "U",
            .l => "L",
            .v => "V",
            .m => "M",
        }, n, alpha, a, lda, px, incx, py, incy),
    }
}

fn mvRun(allocator: Allocator, comptime T: type, comptime name: []const u8, n: BlasInt, alpha: T, a: [*]T, lda: BlasInt, x: [*]T, incx: BlasInt, y: [*]T, incy: BlasInt) Allocator.Error!void {
    const F = Float(T);
    const complex = comptime isComplex(T);
    // symv has no `V` and `M` kernels; `symv` returned above for those
    if (comptime !@hasField(Vtable, letter(T) ++ name)) unreachable;
    // what the core's kernel needs, whichever it is (`buffers.md`)
    const cs = if (complex) 2 else 1;
    const un: usize = @intCast(n);
    const bytes = @max(sym.composeBytes(F, cs, un, incx == 1, incy == 1), sym.sseBytes(F, cs, un));
    // the portable kernel on a core without its own gemv runs plain loops,
    // which need no buffer (`symmetric.Auto`)
    const lower = comptime name[name.len - 1] == 'L' or name[name.len - 1] == 'M';
    const kind: sym.Kind = comptime if (name[0] == 's') .symmetric else if (name[name.len - 1] == 'V' or name[name.len - 1] == 'M') .hermitian_rev else .hermitian;
    const auto = sym.Auto(letter(T), T, lower, kind);
    const direct = @field(blas.table().*, letter(T) ++ name) == (if (complex) &auto.cplx else &auto.real) and !auto.native();
    var none: [1]F = undefined;
    const sc: Scratch(F) = if (direct) .{ .ptr = &none, .len = 0 } else try Scratch(F).heapOnly(allocator, (bytes + @sizeOf(F) - 1) / @sizeOf(F));
    defer sc.deinit(allocator);
    // x86_64's SSE `symv` and `hemv` read past their copies of `x` and `y`,
    // and NaN there reaches `y`
    if (sysinfo.cpu.arch == .x86_64) @memset(sc.slice(), 0);
    // the reference threads big problems here (`symv_thread_*`)
    if (comptime complex) {
        _ = slot(T, name)(long(n), long(n), alpha.re, alpha.im, raw(T, a), long(lda), raw(T, x), incx, raw(T, y), incy, sc.ptr);
    } else {
        _ = slot(T, name)(long(n), long(n), alpha, a, long(lda), x, incx, y, incy, sc.ptr);
    }
}

// ---- tests -------------------------------------------------------------------

const testing = std.testing;
/// Scratch from the heap holds NaN, so a missing clear shows.
var poison: guard.Poison = .{ .child = testing.allocator };
const all = .{ f32, f64, Complex32, Complex64 };

/// Points the global vtable at the portable kernels until
/// `defer blas.testRestore(saved)`.
fn usePortable() blas.TestSaved {
    var vt: Vtable = .{};
    portable.install(&vt);
    return blas.testUse(&vt);
}

fn mulT(comptime T: type, p: T, q: T) T {
    return if (comptime isComplex(T)) p.mul(q) else p * q;
}

fn addT(comptime T: type, p: T, q: T) T {
    return if (comptime isComplex(T)) p.add(q) else p + q;
}

fn conjT(comptime T: type, p: T) T {
    return if (comptime isComplex(T)) p.conjugate() else p;
}

fn zeroT(comptime T: type) T {
    return if (comptime isComplex(T)) .init(0, 0) else 0;
}

fn val(comptime T: type, re: Float(T), im: Float(T)) T {
    return if (comptime isComplex(T)) .init(re, im) else re;
}

/// Element `i` of a vector of `n` at stride `inc`, as an index into its
/// storage: stride -2 puts element 0 last.
fn vi(n: usize, inc: BlasInt, i: usize) usize {
    const k: usize = @abs(inc);
    return if (inc < 0) (n - 1 - i) * k else i * k;
}

fn randomVec(comptime T: type, n: usize, inc: BlasInt, seed: u64) ![]T {
    const v = try testing.allocator.alloc(T, guard.span(n, @abs(inc)));
    check.fill(T, v, seed);
    return v;
}

fn expectSame(comptime T: type, want: []const T, got: []const T, terms: usize) !void {
    try check.expectCloseSlices(T, want, got, check.tolerance(Float(T), 8 * (terms + 2)));
}

const Case = struct { m: usize, n: usize };
const shapes = [_]Case{ .{ .m = 1, .n = 1 }, .{ .m = 5, .n = 3 }, .{ .m = 3, .n = 5 }, .{ .m = 9, .n = 7 }, .{ .m = 17, .n = 4 } };
const stride_pairs = [_][2]BlasInt{ .{ 1, 1 }, .{ 2, -3 }, .{ -1, 2 }, .{ -2, -1 } };

/// `alpha`, `beta` pairs: general, and each special case.
fn scalars(comptime T: type) [6][2]T {
    return .{
        .{ val(T, 1.25, -0.5), val(T, 0.5, 0.25) },
        .{ val(T, 1, 0), val(T, 1, 0) },
        .{ val(T, 0, 0), val(T, 0.5, 0) },
        .{ val(T, 1.5, 0.5), val(T, 0, 0) },
        .{ val(T, 0, 0), val(T, 1, 0) },
        .{ val(T, -1, 0), val(T, 2, 0) },
    };
}

fn testGemv(comptime T: type, trans: Trans, m: usize, n: usize, incx: BlasInt, incy: BlasInt, alpha: T, beta: T) !void {
    const tr = transposed(trans);
    const lenx = if (tr) m else n;
    const leny = if (tr) n else m;
    const lda = m + 2;
    const a = try testing.allocator.alloc(T, lda * n);
    defer testing.allocator.free(a);
    check.fill(T, a, m * 7 + n);
    const x = try randomVec(T, lenx, incx, 3);
    defer testing.allocator.free(x);
    const y = try randomVec(T, leny, incy, 4);
    defer testing.allocator.free(y);
    const want = try testing.allocator.dupe(T, y);
    defer testing.allocator.free(want);
    const conj_a = trans == .r or trans == .c or trans == .s or trans == .d;
    const conj_x = trans == .o or trans == .u or trans == .s or trans == .d;
    for (0..leny) |r| {
        var acc = zeroT(T);
        for (0..lenx) |q| {
            var e = if (tr) a[q + r * lda] else a[r + q * lda];
            var xv = x[vi(lenx, incx, q)];
            if (conj_a) e = conjT(T, e);
            if (conj_x) xv = conjT(T, xv);
            acc = addT(T, acc, mulT(T, e, xv));
        }
        const k = vi(leny, incy, r);
        const scaled = if (isZero(T, beta)) zeroT(T) else mulT(T, beta, y[k]);
        want[k] = if (isZero(T, alpha)) scaled else addT(T, scaled, mulT(T, alpha, acc));
    }
    try gemv(poison.allocator(), T, trans, @intCast(m), @intCast(n), alpha, a.ptr, @intCast(lda), x.ptr, incx, beta, y.ptr, incy);
    try expectSame(T, want, y, lenx);
}

test "interface: gemv matches plain arithmetic, for every trans, stride and scalar case" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (all) |T| {
        const transes = if (comptime isComplex(T)) std.enums.values(Trans) else &[_]Trans{ .n, .t, .r, .c };
        for (transes) |trans| for (shapes) |s| for (stride_pairs) |inc| for (scalars(T)) |ab| {
            testGemv(T, trans, s.m, s.n, inc[0], inc[1], ab[0], ab[1]) catch |e| {
                check.say("gemv {s} {s} m={d} n={d} inc={d},{d}\n", .{ @typeName(T), @tagName(trans), s.m, s.n, inc[0], inc[1] });
                return e;
            };
        };
    }
}

test "interface: gemv beta 0 overwrites NaN in y, and alpha 0 never reads x or a" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (all) |T| {
        const nan = val(T, std.math.nan(Float(T)), std.math.nan(Float(T)));
        var a = [_]T{ nan, nan, nan, nan };
        var x = [_]T{ nan, nan };
        var y = [_]T{ nan, nan };
        try gemv(testing.allocator, T, .n, 2, 2, val(T, 0, 0), &a, 2, &x, 1, val(T, 0, 0), &y, 1);
        try expectSame(T, &[_]T{ zeroT(T), zeroT(T) }, &y, 1);
        y = .{ val(T, 1, 0), val(T, 2, 0) };
        try gemv(testing.allocator, T, .t, 2, 2, val(T, 0, 0), &a, 2, &x, 1, val(T, 2, 0), &y, 1);
        try expectSame(T, &[_]T{ val(T, 2, 0), val(T, 4, 0) }, &y, 1);
    }
}

test "interface: gemv, ger, symv and hemv report the first bad argument and do nothing for empty problems" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    const al = testing.allocator;
    var a = [_]f64{ 1, 2, 3, 4 };
    var x = [_]f64{ 1, 1 };
    var y = [_]f64{ 5, 6 };
    try testing.expectError(error.InvalidTrans, gemv(al, f64, .o, 2, 2, 1, &a, 2, &x, 1, 1, &y, 1));
    try testing.expectError(error.InvalidM, gemv(al, f64, .n, -1, 2, 1, &a, 2, &x, 1, 1, &y, 1));
    try testing.expectError(error.InvalidN, gemv(al, f64, .n, 2, -1, 1, &a, 2, &x, 1, 1, &y, 1));
    try testing.expectError(error.InvalidLda, gemv(al, f64, .n, 2, 2, 1, &a, 1, &x, 1, 1, &y, 1));
    try testing.expectError(error.InvalidIncx, gemv(al, f64, .n, 2, 2, 1, &a, 2, &x, 0, 1, &y, 1));
    try testing.expectError(error.InvalidIncy, gemv(al, f64, .n, 2, 2, 1, &a, 2, &x, 1, 1, &y, 0));
    // m wins over a bad stride
    try testing.expectError(error.InvalidM, gemv(al, f64, .n, -1, 2, 1, &a, 2, &x, 0, 1, &y, 0));
    try testing.expectError(error.InvalidM, ger(al, f64, -1, 2, 1, &x, 1, &y, 1, &a, 2));
    try testing.expectError(error.InvalidN, ger(al, f64, 2, -1, 1, &x, 1, &y, 1, &a, 2));
    try testing.expectError(error.InvalidIncx, ger(al, f64, 2, 2, 1, &x, 0, &y, 1, &a, 2));
    try testing.expectError(error.InvalidIncy, ger(al, f64, 2, 2, 1, &x, 1, &y, 0, &a, 2));
    try testing.expectError(error.InvalidLda, ger(al, f64, 2, 2, 1, &x, 1, &y, 1, &a, 1));
    try testing.expectError(error.InvalidUplo, symv(al, f64, .v, 2, 1, &a, 2, &x, 1, 1, &y, 1));
    try testing.expectError(error.InvalidN, symv(al, f64, .u, -1, 1, &a, 2, &x, 1, 1, &y, 1));
    try testing.expectError(error.InvalidLda, symv(al, f64, .u, 2, 1, &a, 1, &x, 1, 1, &y, 1));
    try testing.expectError(error.InvalidIncx, symv(al, f64, .l, 2, 1, &a, 2, &x, 0, 1, &y, 1));
    try testing.expectError(error.InvalidIncy, symv(al, f64, .l, 2, 1, &a, 2, &x, 1, 1, &y, 0));
    var za = [_]Complex64{ .init(1, 0), .init(0, 0), .init(0, 0), .init(1, 0) };
    var zx = [_]Complex64{ .init(1, 0), .init(1, 0) };
    var zy = [_]Complex64{ .init(1, 0), .init(1, 0) };
    try testing.expectError(error.InvalidLda, hemv(al, Complex64, .m, 2, .init(1, 0), &za, 1, &zx, 1, .init(1, 0), &zy, 1));
    // empty problems leave y alone, even with a scaling beta
    try gemv(al, f64, .n, 0, 2, 1, &a, 1, &x, 1, 2, &y, 1);
    try gemv(al, f64, .t, 2, 0, 1, &a, 2, &x, 1, 2, &y, 1);
    try ger(al, f64, 0, 2, 1, &x, 1, &y, 1, &a, 1);
    try symv(al, f64, .u, 0, 1, &a, 1, &x, 1, 2, &y, 1);
    try testing.expectEqualSlices(f64, &.{ 5, 6 }, &y);
    try testing.expectEqualSlices(f64, &.{ 1, 2, 3, 4 }, &a);
}

fn testGer(comptime T: type, conj: bool, m: usize, n: usize, incx: BlasInt, incy: BlasInt, alpha: T) !void {
    const lda = m + 2;
    const a = try testing.allocator.alloc(T, lda * n);
    defer testing.allocator.free(a);
    check.fill(T, a, m * 5 + n);
    const x = try randomVec(T, m, incx, 6);
    defer testing.allocator.free(x);
    const y = try randomVec(T, n, incy, 7);
    defer testing.allocator.free(y);
    const want = try testing.allocator.dupe(T, a);
    defer testing.allocator.free(want);
    for (0..n) |j| for (0..m) |i| {
        const yv = y[vi(n, incy, j)];
        const t = mulT(T, x[vi(m, incx, i)], if (conj) conjT(T, yv) else yv);
        want[i + j * lda] = addT(T, want[i + j * lda], mulT(T, alpha, t));
    };
    if (conj) {
        if (comptime !isComplex(T)) unreachable;
        try gerc(poison.allocator(), T, @intCast(m), @intCast(n), alpha, x.ptr, incx, y.ptr, incy, a.ptr, @intCast(lda));
    } else {
        try ger(poison.allocator(), T, @intCast(m), @intCast(n), alpha, x.ptr, incx, y.ptr, incy, a.ptr, @intCast(lda));
    }
    try expectSame(T, want, a, 1);
}

test "interface: ger and gerc match plain arithmetic, for every stride and an alpha of 0" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (all) |T| {
        for (shapes) |s| for (stride_pairs) |inc| for ([_]T{ val(T, 1.25, -0.5), val(T, 0, 0), val(T, -1, 0) }) |alpha| {
            testGer(T, false, s.m, s.n, inc[0], inc[1], alpha) catch |e| {
                check.say("ger {s} m={d} n={d} inc={d},{d}\n", .{ @typeName(T), s.m, s.n, inc[0], inc[1] });
                return e;
            };
            if (comptime isComplex(T)) try testGer(T, true, s.m, s.n, inc[0], inc[1], alpha);
        };
    }
}

/// `S[i][j]` of the symmetric or Hermitian matrix `uplo` stores in `a`.
fn symAt(comptime T: type, kind: sym.Kind, upper: bool, a: []const T, lda: usize, i: usize, j: usize) T {
    const stored = if (upper) i <= j else i >= j;
    const e = if (stored) a[i + j * lda] else a[j + i * lda];
    if (i == j) return if (comptime !isComplex(T)) e else if (kind == .symmetric) e else val(T, e.re, 0);
    return switch (kind) {
        .symmetric => e,
        .hermitian => if (stored) e else conjT(T, e),
        .hermitian_rev => if (stored) conjT(T, e) else e,
    };
}

fn testMv(comptime T: type, herm: bool, uplo: Uplo, n: usize, incx: BlasInt, incy: BlasInt, alpha: T, beta: T) !void {
    const lda = n + 1;
    const a = try testing.allocator.alloc(T, lda * n);
    defer testing.allocator.free(a);
    check.fill(T, a, n * 3 + 1);
    const x = try randomVec(T, n, incx, 8);
    defer testing.allocator.free(x);
    const y = try randomVec(T, n, incy, 9);
    defer testing.allocator.free(y);
    const want = try testing.allocator.dupe(T, y);
    defer testing.allocator.free(want);
    const upper = uplo == .u or uplo == .v;
    const kind: sym.Kind = if (uplo == .v or uplo == .m) .hermitian_rev else if (herm) .hermitian else .symmetric;
    for (0..n) |i| {
        var acc = zeroT(T);
        for (0..n) |j| acc = addT(T, acc, mulT(T, symAt(T, kind, upper, a, lda, i, j), x[vi(n, incx, j)]));
        const k = vi(n, incy, i);
        const scaled = if (isZero(T, beta)) zeroT(T) else mulT(T, beta, y[k]);
        want[k] = if (isZero(T, alpha)) scaled else addT(T, scaled, mulT(T, alpha, acc));
    }
    if (kind == .symmetric) {
        try symv(poison.allocator(), T, uplo, @intCast(n), alpha, a.ptr, @intCast(lda), x.ptr, incx, beta, y.ptr, incy);
    } else if (comptime isComplex(T)) {
        try hemv(poison.allocator(), T, uplo, @intCast(n), alpha, a.ptr, @intCast(lda), x.ptr, incx, beta, y.ptr, incy);
    }
    try expectSame(T, want, y, n);
}

const mv_sizes = [_]usize{ 1, 2, 5, 8, 9, 17 };

test "interface: symv and hemv match plain arithmetic, for every triangle, stride and scalar case" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (all) |T| {
        for ([_]bool{ false, true }) |herm| {
            if (herm and !comptime isComplex(T)) continue;
            const uplos = if (herm) std.enums.values(Uplo) else &[_]Uplo{ .u, .l };
            for (uplos) |uplo| for (mv_sizes) |n| for (stride_pairs) |inc| for (scalars(T)) |ab| {
                testMv(T, herm, uplo, n, inc[0], inc[1], ab[0], ab[1]) catch |e| {
                    check.say("{s} {s} {s} n={d} inc={d},{d}\n", .{ if (herm) "hemv" else "symv", @typeName(T), @tagName(uplo), n, inc[0], inc[1] });
                    return e;
                };
            };
        }
    }
}

test "interface: problems over the stack buffer allocate, and small gemv and ger never do" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (all) |T| {
        try testGemv(T, .n, 300, 300, 1, 1, val(T, 1.5, 0.5), val(T, 0.5, 0));
        try testGemv(T, if (isComplex(T)) .c else .t, 300, 300, -2, 3, val(T, 1.5, 0.5), val(T, 1, 0));
        try testGer(T, false, 600, 3, 2, -1, val(T, 0.75, 0.25));
        try testMv(T, false, .l, 300, -1, 2, val(T, 1.5, 0.5), val(T, 0.5, 0));
        if (comptime isComplex(T)) {
            try testMv(T, true, .u, 300, 2, 1, val(T, 1.5, 0.5), val(T, 0.5, 0));
        }
        // a failing allocator is never asked for tiny problems
        var a = [_]T{ val(T, 1, 0), val(T, 2, 0), val(T, 3, 0), val(T, 4, 0) };
        var x = [_]T{ val(T, 1, 0), val(T, 1, 0) };
        var y = [_]T{ val(T, 1, 0), val(T, 1, 0) };
        const fa = testing.failing_allocator;
        try gemv(fa, T, .n, 2, 2, val(T, 1, 0), &a, 2, &x, 1, val(T, 1, 0), &y, 1);
        try ger(fa, T, 2, 2, val(T, 1, 0), &x, 2, &y, 1, &a, 2);
        // nor for symv on the portable gemv (plain loops), or a unit-stride ger
        try symv(fa, T, .u, 2, val(T, 1, 0), &a, 2, &x, 1, val(T, 1, 0), &y, 1);
        const big = 600;
        const bx = try testing.allocator.alloc(T, 2 * big);
        defer testing.allocator.free(bx);
        check.fill(T, bx, 3);
        try ger(fa, T, big, 1, val(T, 1, 0), bx.ptr, 1, &y, 2, bx.ptr, big);
        // and a big one that needs a buffer reports the failure
        try testing.expectError(error.OutOfMemory, ger(fa, T, big, 1, val(T, 1, 0), bx.ptr, 2, &y, 2, bx.ptr, big));
    }
}

// ---- half-precision gemv tests ----------------------------------------------

const bf = @import("../ext/bf16/__root__.zig");

fn halfIn(comptime E: type, v: f32) E {
    return if (E == types.BF16) bf.toBf16(v) else @floatCast(v);
}

fn halfVal(comptime E: type, v: E) f64 {
    return if (E == types.BF16) bf.fromBf16(v) else @floatCast(v);
}

fn testGemvHalf(comptime kind: HalfGemv, trans: Trans, m: usize, n: usize, incx: BlasInt, incy: BlasInt, alpha: f32, beta: f32) !void {
    const In = kind.In();
    const Out = kind.Out();
    const tr = transposed(trans);
    const lenx = if (tr) m else n;
    const leny = if (tr) n else m;
    const lda = m + 2;
    const a = try testing.allocator.alloc(In, lda * n);
    defer testing.allocator.free(a);
    const x = try testing.allocator.alloc(In, guard.span(lenx, @abs(incx)));
    defer testing.allocator.free(x);
    const y = try testing.allocator.alloc(Out, guard.span(leny, @abs(incy)));
    defer testing.allocator.free(y);
    var f: [1024]f32 = undefined;
    check.fill(f32, f[0..a.len], m * 7 + n);
    for (a, f[0..a.len]) |*d, v| d.* = halfIn(In, v);
    check.fill(f32, f[0..x.len], 3);
    for (x, f[0..x.len]) |*d, v| d.* = halfIn(In, v);
    check.fill(f32, f[0..y.len], 4);
    for (y, f[0..y.len]) |*d, v| d.* = halfIn(Out, v);
    const want = try testing.allocator.dupe(Out, y);
    defer testing.allocator.free(want);
    for (0..leny) |r| {
        var acc: f64 = 0;
        for (0..lenx) |q| {
            const e = if (tr) a[q + r * lda] else a[r + q * lda];
            acc += halfVal(In, e) * halfVal(In, x[vi(lenx, incx, q)]);
        }
        const k = vi(leny, incy, r);
        const af: f64 = halfVal(Out, halfIn(Out, alpha));
        const bt: f64 = halfVal(Out, halfIn(Out, beta));
        const prior: f64 = if (bt == 0) 0 else bt * halfVal(Out, y[k]);
        want[k] = halfIn(Out, @floatCast(if (af == 0) prior else af * acc + prior));
    }
    try gemvHalf(kind, trans, @intCast(m), @intCast(n), halfIn(Out, alpha), a.ptr, @intCast(lda), x.ptr, incx, halfIn(Out, beta), y.ptr, incy);
    for (want, y, 0..) |w, g, i| {
        const wv = halfVal(Out, w);
        const tol: f64 = if (Out == types.BF16) 0.012 * @abs(wv) + 1e-3 else 1e-5 * @as(f64, @floatFromInt(lenx + 2)) + 1e-6 * @abs(wv);
        if (!(@abs(wv - halfVal(Out, g)) <= tol)) {
            check.say("{s} {s} m={d} n={d} inc={d},{d} alpha={d} beta={d} y[{d}]: got {d}, want {d}\n", .{ @tagName(kind), @tagName(trans), m, n, incx, incy, alpha, beta, i, halfVal(Out, g), wv });
            return error.TestExpectedApproxEqAbs;
        }
    }
}

const half_scalars = [_][2]f32{ .{ 1.25, 0.5 }, .{ 1, 1 }, .{ 0, 0.5 }, .{ 1.5, 0 }, .{ 0, 1 }, .{ -1, 2 }, .{ 0, 0 } };

test "interface: the half-precision gemvs match plain arithmetic, for every trans, stride and scalar case" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (std.enums.values(HalfGemv)) |kind| {
        for ([_]Trans{ .n, .t, .r, .c }) |trans| for (shapes) |s| for (stride_pairs) |inc| for (half_scalars) |ab| {
            testGemvHalf(kind, trans, s.m, s.n, inc[0], inc[1], ab[0], ab[1]) catch |e| return e;
        };
    }
}

test "interface: half-precision gemv beta 0 overwrites NaN, alpha 0 never reads x or a" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (std.enums.values(HalfGemv)) |kind| {
        const In = kind.In();
        const Out = kind.Out();
        const nan_in = halfIn(In, std.math.nan(f32));
        var a = [_]In{ nan_in, nan_in, nan_in, nan_in };
        var x = [_]In{ nan_in, nan_in };
        var y = [_]Out{ halfIn(Out, std.math.nan(f32)), halfIn(Out, std.math.nan(f32)) };
        try gemvHalf(kind, .n, 2, 2, halfIn(Out, 0), &a, 2, &x, 1, halfIn(Out, 0), &y, 1);
        try testing.expectEqual(@as(f64, 0), halfVal(Out, y[0]));
        try testing.expectEqual(@as(f64, 0), halfVal(Out, y[1]));
        y = .{ halfIn(Out, 1), halfIn(Out, 2) };
        try gemvHalf(kind, .t, 2, 2, halfIn(Out, 0), &a, 2, &x, -1, halfIn(Out, 2), &y, 1);
        try testing.expectEqual(@as(f64, 2), halfVal(Out, y[0]));
        try testing.expectEqual(@as(f64, 4), halfVal(Out, y[1]));
        // beta 1 with alpha 0 touches nothing, even through a stride
        y = .{ halfIn(Out, 1), halfIn(Out, 2) };
        try gemvHalf(kind, .n, 2, 2, halfIn(Out, 0), &a, 2, &x, 1, halfIn(Out, 1), &y, -1);
        try testing.expectEqual(@as(f64, 1), halfVal(Out, y[0]));
    }
}

test "interface: half-precision gemv reports the first bad argument and ignores empty problems" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (std.enums.values(HalfGemv)) |kind| {
        const In = kind.In();
        const Out = kind.Out();
        var a = [_]In{ halfIn(In, 1), halfIn(In, 2), halfIn(In, 3), halfIn(In, 4) };
        var x = [_]In{ halfIn(In, 1), halfIn(In, 1) };
        var y = [_]Out{ halfIn(Out, 5), halfIn(Out, 6) };
        const one = halfIn(Out, 1);
        try testing.expectError(error.InvalidTrans, gemvHalf(kind, .o, 2, 2, one, &a, 2, &x, 1, one, &y, 1));
        try testing.expectError(error.InvalidTrans, gemvHalf(kind, .d, -1, -1, one, &a, 0, &x, 0, one, &y, 0));
        try testing.expectError(error.InvalidM, gemvHalf(kind, .n, -1, -1, one, &a, 0, &x, 0, one, &y, 0));
        try testing.expectError(error.InvalidN, gemvHalf(kind, .n, 2, -1, one, &a, 0, &x, 0, one, &y, 0));
        try testing.expectError(error.InvalidLda, gemvHalf(kind, .n, 2, 2, one, &a, 1, &x, 0, one, &y, 0));
        try testing.expectError(error.InvalidIncx, gemvHalf(kind, .t, 2, 2, one, &a, 2, &x, 0, one, &y, 0));
        try testing.expectError(error.InvalidIncy, gemvHalf(kind, .t, 2, 2, one, &a, 2, &x, 1, one, &y, 0));
        // an empty problem does not even apply beta
        try gemvHalf(kind, .n, 0, 2, one, &a, 1, &x, 1, halfIn(Out, 0), &y, 1);
        try gemvHalf(kind, .n, 2, 0, one, &a, 2, &x, 1, halfIn(Out, 0), &y, 1);
        try testing.expectEqual(@as(f64, 5), halfVal(Out, y[0]));
    }
}

fn checkAll(vt: *const Vtable) !void {
    // the interface reaches kernels through the global vtable
    const saved = blas.testUse(vt);
    defer blas.testRestore(saved);
    inline for (std.enums.values(HalfGemv)) |kind| {
        try testGemvHalf(kind, .n, 9, 7, 2, 1, 1.25, 0.5);
        try testGemvHalf(kind, .t, 9, 7, -1, 3, 1.25, 0.5);
        try testGemvHalf(kind, .n, 33, 4, 1, 1, 1.25, 0);
        try testGemvHalf(kind, .t, 4, 33, 1, 1, 1.25, 1);
        try testGemvHalf(kind, .r, 17, 3, 1, -2, 0, 0.5);
    }
    inline for (all) |T| {
        try testGemv(T, .n, 9, 7, 2, 1, val(T, 1.25, -0.5), val(T, 0.5, 0.25));
        try testGemv(T, .t, 9, 7, -1, 3, val(T, 1.25, -0.5), val(T, 0.5, 0.25));
        try testGemv(T, .n, 33, 4, 1, 1, val(T, 1.25, -0.5), val(T, 0, 0));
        try testGemv(T, .t, 4, 33, 1, 1, val(T, 1.25, -0.5), val(T, 1, 0));
        if (comptime isComplex(T)) {
            inline for (.{ Trans.r, Trans.c, Trans.o, Trans.u, Trans.s, Trans.d }) |t| try testGemv(T, t, 9, 7, 2, -1, val(T, 1.25, -0.5), val(T, 0.5, 0.25));
        }
        try testGer(T, false, 17, 9, 1, 1, val(T, 1.25, -0.5));
        try testGer(T, false, 9, 17, -2, 3, val(T, 1.25, -0.5));
        if (comptime isComplex(T)) try testGer(T, true, 9, 17, 2, 1, val(T, 1.25, -0.5));
        for ([_]bool{ false, true }) |herm| {
            if (herm and !comptime isComplex(T)) continue;
            for (if (herm) std.enums.values(Uplo) else &[_]Uplo{ .u, .l }) |uplo| {
                for ([_]usize{ 1, 7, 8, 17, 33 }) |n| {
                    try testMv(T, herm, uplo, n, 1, 1, val(T, 1.25, -0.5), val(T, 0.5, 0.25));
                    try testMv(T, herm, uplo, n, -2, 3, val(T, 1.25, -0.5), val(T, 1, 0));
                }
            }
        }
    }
}

test "every linked core's kernels work through the level-2 interface, where this machine runs them" {
    try check.eachCore(checkAll);
}
