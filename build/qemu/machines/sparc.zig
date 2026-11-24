//! The SPARC rigs.

const std = @import("std");
const machines = @import("../machines.zig");

const Machine = machines.Machine;
const small_filters = machines.small_filters;
const sparc64_stub = machines.sparc64_stub;

pub const list = [_]Machine{
    .{
        .name = "sparc64",
        .stub = sparc64_stub,
        .optimize = .fast,
        // OpenBIOS maps 0x400000 to 0x800000 only, and the default set's
        // code alone is about 5 MB; blas2 by itself is 3.3 MB
        .parts = &.{
            .{ .name = "blas1", .filters = &.{ "math.kernels.blas.blas1.", "math.kernels.blas._test.blas1.", "math.kernels.blas._test.check.", "math.kernels.blas._test.bitwise." } },
            .{ .name = "blas2-dense", .filters = &.{ "math.kernels.blas.blas2.__root__.", "math.kernels.blas.blas2.gemv.", "math.kernels.blas._test.blas2.gemv.", "math.kernels.blas._test.bf16.", "math.kernels.blas.blas2.ger.", "math.kernels.blas.blas2.hemv.", "math.kernels.blas.blas2.symv." } },
            .{ .name = "blas2-banded", .filters = &.{ "math.kernels.blas._test.blas2.banded.", "math.kernels.blas.blas2.symmetric." } },
            .{ .name = "blas2-triangular", .filters = &.{"math.kernels.blas._test.blas2.triangular."} },
            .{ .name = "blas2-interface", .filters = &.{ "math.kernels.blas._test.blas2.interface.", "math.kernels.blas._test.blas2.interface_any." } },
            .{ .name = "blas2-edge", .filters = &.{"math.kernels.blas._test.blas2.edge."} },
            // one per routine family: together past the window
            .{ .name = "blas3-gemm", .filters = &.{ "math.kernels.blas.blas3.__root__.", "math.kernels.blas.blas3.gemm.", "math.kernels.blas._test.blas3.params." } },
            .{ .name = "blas3-symm", .filters = &.{"math.kernels.blas.blas3.symm."} },
            .{ .name = "blas3-rank", .filters = &.{ "math.kernels.blas.blas3.syrk.", "math.kernels.blas.blas3.syr2k." } },
            .{ .name = "blas3-tri", .filters = &.{ "math.kernels.blas.blas3.trmm.", "math.kernels.blas.blas3.trsm.", "math.kernels.lapack." } },
            .{ .name = "blas3-interface", .filters = &.{"math.kernels.blas._test.blas3.interface."} },
            .{ .name = "openblas", .filters = &.{"math.kernels.blas._test.from_openblas."} },
            .{ .name = "rest", .filters = &.{ "math.kernels.blas.__root__.", "math.kernels.blas.Config.", "math.kernels.blas.cpu." } },
            .{ .name = "small", .filters = &small_filters },
        },
        .query = .{ .cpu_arch = .sparc64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.sparc.cpu.v9 } },
        .qemu = "qemu-system-sparc64",
        .args = &.{ "-M", "sun4u" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/sparc64.ld",
            .load_window = .{ .base = 0x400000, .size = 0x400000 },
        } },
    },
    .{
        // sparc64 in debug, the mode the other rigs build in; its images are
        // larger, so blas3 is split too, and gemm's three test files apart
        .name = "sparc64-debug",
        .stub = sparc64_stub,
        .optimize = .debug,
        .parts = &.{
            // blas1 in three, as two are past the window
            .{ .name = "blas1", .filters = &.{ "math.kernels.blas.blas1.__root__.", "math.kernels.blas.blas1.asum.", "math.kernels.blas.blas1.axpby.", "math.kernels.blas.blas1.axpy.", "math.kernels.blas._test.check.", "math.kernels.blas._test.bitwise.", "math.kernels.blas.blas1.copy.", "math.kernels.blas.blas1.dot.", "math.kernels.blas.blas1.extremum.", "math.kernels.blas._test.blas1.extremum.", "math.kernels.blas.blas1.nrm2.", "math.kernels.blas.blas1.rot.", "math.kernels.blas.blas1.rotm.", "math.kernels.blas.blas1.scal.", "math.kernels.blas.blas1.sum.", "math.kernels.blas.blas1.swap.", "math.kernels.blas.uk_block." } },
            .{ .name = "blas1-interface", .filters = &.{ "math.kernels.blas._test.blas1.interface.", "math.kernels.blas._test.blas1.interface_any.", "math.kernels.blas._test.blas1.rot_reals.", "math.kernels.blas._test.blas1.threaded." } },
            .{ .name = "blas1-edge", .filters = &.{ "math.kernels.blas._test.blas1.edge.test.portable level-1 kernels: special", "math.kernels.blas._test.blas1.edge.test.every linked core's level-1 kernels: special", "math.kernels.blas._test.blas1.edge.test.portable level-1 reductions", "math.kernels.blas._test.blas1.edge.test.every linked core's level-1 reductions", "math.kernels.blas._test.blas1.edge.test.every linked core keeps" } },
            .{ .name = "blas1-denormal", .filters = &.{ "math.kernels.blas._test.blas1.edge.test.portable level-1 kernels: denormal", "math.kernels.blas._test.blas1.edge.test.every linked core's level-1 kernels: denormal", "math.kernels.blas._test.blas1.edge.test.where util.vectorIeee" } },
            // the dense routines, the triangular ones, the interface and the
            // complex edge cases each in two, as one is past the window
            .{ .name = "blas2-dense", .filters = &.{ "math.kernels.blas.blas2.__root__.", "math.kernels.blas._test.blas2.gemv.test.portable", "math.kernels.blas.blas2.gemv.__root__.test_" } },
            .{ .name = "blas2-gemv-x87", .filters = &.{ "math.kernels.blas._test.bf16.gemv.", "math.kernels.blas.blas2.gemv.extended.test.extended gemv:" } },
            .{ .name = "blas2-gemv-cores", .filters = &.{"math.kernels.blas._test.blas2.gemv.test.every linked core's"} },
            .{ .name = "blas2-dense-2", .filters = &.{ "math.kernels.blas.blas2.ger.", "math.kernels.blas.blas2.hemv.", "math.kernels.blas.blas2.symv." } },
            .{ .name = "blas2-banded", .filters = &.{
                "math.kernels.blas._test.blas2.banded.test.banded",
                "math.kernels.blas._test.blas2.banded.test.the small",
                "math.kernels.blas._test.blas2.banded.test.strided",
                "math.kernels.blas._test.blas2.banded.test.an allocation",
                "math.kernels.blas._test.blas2.banded.test.invalid",
                "math.kernels.blas.blas2.symmetric.",
            } },
            .{ .name = "blas2-triangular", .filters = &.{"math.kernels.blas._test.blas2.triangular.test.triangular"} },
            // every core's triangular results in three, as one is past the window
            .{ .name = "blas2-triangular-cores-real", .filters = &.{"math.kernels.blas._test.blas2.triangular.test.every linked core's kernels give the same triangular results, real"} },
            .{ .name = "blas2-triangular-cores-complex", .filters = &.{"math.kernels.blas._test.blas2.triangular.test.every linked core's kernels give the same triangular results, complex products"} },
            .{ .name = "blas2-triangular-cores-complex-2", .filters = &.{"math.kernels.blas._test.blas2.triangular.test.every linked core's kernels give the same triangular results, complex solves"} },
            .{ .name = "blas2-f80", .filters = &.{ "math.kernels.blas._test.blas2.banded.test.f80", "math.kernels.blas._test.blas2.triangular.test.f80" } },
            .{ .name = "blas2-interface", .filters = &.{
                "math.kernels.blas._test.blas2.interface.test.interface: gemv",
                "math.kernels.blas._test.blas2.interface.test.interface: ger",
                "math.kernels.blas._test.blas2.interface.test.interface: symv",
                "math.kernels.blas._test.blas2.interface.test.interface: problems",
            } },
            .{ .name = "blas2-interface-2", .filters = &.{
                "math.kernels.blas._test.blas2.interface.test.interface: the half",
                "math.kernels.blas._test.blas2.interface.test.interface: half",
                "math.kernels.blas._test.blas2.interface_any.test.every linked core's",
                "math.kernels.blas._test.blas2.interface.test_",
            } },
            // f16, f80 and f128: past the window with the rest
            .{ .name = "blas2-interface-any", .filters = &.{"math.kernels.blas._test.blas2.interface_any.test.interface: level 2"} },
            // real, complex, and complex on every core in two, as one is past
            // the window
            .{ .name = "blas2-edge-real", .filters = &.{ "math.kernels.blas._test.blas2.edge.test.portable level-2 routines, real", "math.kernels.blas._test.blas2.edge.test.every linked core's level-2 routines, real" } },
            .{ .name = "blas2-edge-complex", .filters = &.{"math.kernels.blas._test.blas2.edge.test.portable level-2 routines, complex"} },
            .{ .name = "blas2-edge-complex-cores", .filters = &.{"math.kernels.blas._test.blas2.edge.test.every linked core's level-2 routines, complex products"} },
            .{ .name = "blas2-edge-complex-cores-2", .filters = &.{"math.kernels.blas._test.blas2.edge.test.every linked core's level-2 routines, complex solves"} },
            .{ .name = "blas3-gemm", .filters = &.{"math.kernels.blas.blas3.gemm.__root__."} },
            // but "every core of" each family: those are the blas3-params-* parts
            .{ .name = "blas3-params", .filters = &.{
                "math.kernels.blas._test.blas3.params.test.dtb entries",
                "math.kernels.blas._test.blas3.params.test.generic follows",
                "math.kernels.blas._test.blas3.params.test.the block factor",
                "math.kernels.blas._test.blas3.params.test.fixed sets",
                "math.kernels.blas._test.blas3.params.test.old x86",
                "math.kernels.blas._test.blas3.params.test.LA464",
                "math.kernels.blas._test.blas3.params.test.detection gives",
            } },
            .{ .name = "blas3-params-x86_64", .filters = &.{"math.kernels.blas._test.blas3.params.test.every core of the x86_64 family"} },
            .{ .name = "blas3-params-x86", .filters = &.{"math.kernels.blas._test.blas3.params.test.every core of the x86 family"} },
            .{ .name = "blas3-params-arm", .filters = &.{ "math.kernels.blas._test.blas3.params.test.every core of the arm64 family", "math.kernels.blas._test.blas3.params.test.every core of the arm family" } },
            .{ .name = "blas3-params-power", .filters = &.{ "math.kernels.blas._test.blas3.params.test.every core of the power64 family", "math.kernels.blas._test.blas3.params.test.every core of the power32 family" } },
            .{ .name = "blas3-params-rest", .filters = &.{ "math.kernels.blas._test.blas3.params.test.every core of the riscv64 family", "math.kernels.blas._test.blas3.params.test.every core of the s390x family", "math.kernels.blas._test.blas3.params.test.every core of the mips64 family", "math.kernels.blas._test.blas3.params.test.every core of the mips32 family", "math.kernels.blas._test.blas3.params.test.every core of the loongarch64 family", "math.kernels.blas._test.blas3.params.test.every core of the sparc family", "math.kernels.blas._test.blas3.params.test.every core of the wasm family", "math.kernels.blas._test.blas3.params.test.every core of the generic family" } },
            .{ .name = "blas3-kernel", .filters = &.{
                "math.kernels.blas.blas3.__root__.",
                "math.kernels.blas.blas3.gemm.kernel.test.packed kernels",
                "math.kernels.blas.blas3.gemm.kernel.test.geometry",
                "math.kernels.blas.blas3.gemm.kernel.test.beta",
            } },
            // the complex kernels in three, as the 8-lane tile's alone are
            // past the window
            .{ .name = "blas3-kernel-complex", .filters = &.{"math.kernels.blas.blas3.gemm.kernel.test.complex packed kernels match"} },
            .{ .name = "blas3-kernel-complex-wide", .filters = &.{"math.kernels.blas.blas3.gemm.kernel.test.complex packed kernels of an 8-lane tile match plain arithmetic, A plain"} },
            .{ .name = "blas3-kernel-complex-wide-2", .filters = &.{"math.kernels.blas.blas3.gemm.kernel.test.complex packed kernels of an 8-lane tile match plain arithmetic, A conjugated"} },
            .{ .name = "blas3-symm", .filters = &.{"math.kernels.blas.blas3.symm."} },
            .{ .name = "blas3-syrk", .filters = &.{"math.kernels.blas.blas3.syrk."} },
            .{ .name = "blas3-syr2k", .filters = &.{"math.kernels.blas.blas3.syr2k."} },
            .{ .name = "blas3-trmm", .filters = &.{"math.kernels.blas.blas3.trmm."} },
            .{ .name = "blas3-trsm", .filters = &.{"math.kernels.blas.blas3.trsm."} },
            // the interface in five, as two are near the window
            .{ .name = "blas3-interface", .filters = &.{ "math.kernels.blas._test.blas3.interface.test.interface: gemm", "math.kernels.blas._test.blas3.interface.test.interface: a matrix past" } },
            .{ .name = "blas3-interface-symm", .filters = &.{"math.kernels.blas._test.blas3.interface.test.interface: symm"} },
            .{ .name = "blas3-interface-rank", .filters = &.{ "math.kernels.blas._test.blas3.interface.test.interface: syrk", "math.kernels.blas._test.blas3.interface.test.interface: the rank-k" } },
            .{ .name = "blas3-interface-tri", .filters = &.{"math.kernels.blas._test.blas3.interface.test.interface: trmm"} },
            .{ .name = "blas3-interface-tri-unslotted", .filters = &.{"math.kernels.blas._test.blas3.interface.test.interface: unslotted trmm"} },
            // the OpenBLAS suites in three
            .{ .name = "openblas-l1", .filters = &.{ "math.kernels.blas._test.from_openblas.blat1.", "math.kernels.blas._test.from_openblas.utest.", "math.kernels.blas._test.from_openblas.utest_ext." } },
            // blat2 in two, past the window with level 3's kernels
            .{ .name = "openblas-blat2", .filters = &.{ "math.kernels.blas._test.from_openblas.blat2.test_", "math.kernels.blas._test.from_openblas.blat2.test.blat2: SMVCH", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xGEMV", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xGBMV", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xSYMV", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xSBMV", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xSPMV", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xTRMV", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xTBMV", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xTPMV", "math.kernels.blas._test.from_openblas.blat2.test.blat2: error exits" } },
            .{ .name = "openblas-blat2-b", .filters = &.{ "math.kernels.blas._test.from_openblas.blat2.test.blat2: xTRSV", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xTBSV", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xTPSV", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xGER,", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xGERC", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xSYR,", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xSPR,", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xSYR2", "math.kernels.blas._test.from_openblas.blat2.test.blat2: xSPR2" } },
            .{ .name = "openblas-blat3", .filters = &.{"math.kernels.blas._test.from_openblas.blat3."} },
            // the suites on every linked core, by level, level 1 by suite
            // (utest_ext's level 3 runs with level 3), level 2 by type
            .{ .name = "openblas-cores-l1-blat1", .filters = &.{"math.kernels.blas._test.from_openblas.__root__.test.blat1 on every"} },
            .{ .name = "openblas-cores-l1-utest", .filters = &.{"math.kernels.blas._test.from_openblas.__root__.test.utest and utest_ext's level-1 cases"} },
            .{ .name = "openblas-cores-l1-utest-2", .filters = &.{"math.kernels.blas._test.from_openblas.__root__.test.utest_ext's level-2 cases"} },
            .{ .name = "openblas-cores-l2-real", .filters = &.{ "math.kernels.blas._test.from_openblas.__root__.test.blat2 at f32", "math.kernels.blas._test.from_openblas.__root__.test.blat2 at f64" } },
            .{ .name = "openblas-cores-l2-complex", .filters = &.{"math.kernels.blas._test.from_openblas.__root__.test.blat2 at complex"} },
            .{ .name = "openblas-cores-l3", .filters = &.{"math.kernels.blas._test.from_openblas.__root__.test.blat3, utest's LAPACK cases"} },
            // LAPACK in five, as one is past the window
            .{ .name = "lapack-cholesky", .filters = &.{"math.kernels.lapack._test.cholesky."} },
            .{ .name = "lapack-getrf", .filters = &.{ "math.kernels.lapack._test.lu.test.getrf", "math.kernels.lapack._test.lu_solve.test.getrf" } },
            .{ .name = "lapack-getrs-real", .filters = &.{ "math.kernels.lapack._test.lu_solve.test.getrs: op(A) X = B on getrf's factor, every op, NaN in one column of B, real", "math.kernels.lapack._test.lu_solve.test.getrs: a known", "math.kernels.lapack._test.lu_solve.test.getrs: reports" } },
            .{ .name = "lapack-getrs-complex", .filters = &.{"math.kernels.lapack._test.lu_solve.test.getrs: op(A) X = B on getrf's factor, every op, NaN in one column of B, complex"} },
            .{ .name = "lapack-faer", .filters = &.{ "math.kernels.lapack.ops.", "math.kernels.lapack.interface.", "math.kernels.lapack.blas_calls.", "math.kernels.lapack.auxiliary.", "math.kernels.lapack.householder.", "math.kernels.lapack._test.qr.", "math.kernels.lapack.qr_solve.", "math.kernels.lapack.tridiag.", "math.kernels.lapack._test.tridiag_evd.", "math.kernels.lapack.hessenberg.", "math.kernels.lapack.real_schur.", "math.kernels.lapack.complex_schur.", "math.kernels.lapack._test.balance.", "math.kernels.lapack.bidiag.", "math.kernels.lapack._test.bidiag_svd." } },
            .{ .name = "rest", .filters = &.{ "math.kernels.blas.__root__.", "math.kernels.blas.Config.", "math.kernels.blas.cpu." } },
            .{ .name = "small", .filters = &small_filters },
        },
        .query = .{ .cpu_arch = .sparc64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.sparc.cpu.v9 } },
        .qemu = "qemu-system-sparc64",
        .args = &.{ "-M", "sun4u" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/sparc64.ld",
            .load_window = .{ .base = 0x400000, .size = 0x400000 },
        } },
    },
    .{
        // OpenBIOS on sun4m, as sparc64's on sun4u
        .name = "sparc",
        .query = .{ .cpu_arch = .sparc, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.sparc.cpu.v8 } },
        .qemu = "qemu-system-sparc",
        .args = &.{ "-M", "SS-5" },
        .memory = "256M",
        .stub = "lld cannot link 32-bit SPARC, and Zig's own ELF linker crashes on it (issue #189)",
        .boot = .{ .kernel = .{
            .linker_script = "",
        } },
    },
};
