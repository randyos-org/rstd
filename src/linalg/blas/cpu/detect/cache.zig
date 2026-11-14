//! The cache sizes the gemm blocking rules read (`blas3/gemm/params.zig`),
//! where an instruction gives them: cpuid on x86 (`get_l2_size` /
//! `get_l3_size`) and cpucfg on LoongArch. Elsewhere, and where the answer
//! is zero, the defaults stand (OpenBLAS assumes 256 KB on x86 and 512 KB on
//! RISC-V, which reads sysfs and so needs an OS).

const sysinfo = @import("builtin");
const x86 = @import("x86.zig");
const loongarch = @import("loongarch.zig");

pub const Cache = struct {
    /// KB.
    l2_kb: i32 = 512,
    l3_kb: i32 = 0,
    /// x86 only: an AMD part with AVX-512 (Zen 4 and later).
    amd_avx512: bool = false,
    /// LoongArch processor id (`cpucfg` 0).
    prid: u32 = 0,
};

pub fn detect() Cache {
    var c: Cache = .{};
    switch (sysinfo.target.cpu.arch) {
        .x86_64, .x86 => {
            const id = x86.cpuid(0x80000006, 0);
            const l2: i32 = @intCast((id.ecx >> 16) & 0xffff);
            c.l2_kb = if (l2 > 0) l2 else 256;
            c.l3_kb = @intCast(((id.edx >> 18) & 0x3fff) * 512);
            const vendor = x86.cpuid(0, 0);
            const leaf7 = x86.cpuid(7, 0);
            // "Auth" of "AuthenticAMD"
            c.amd_avx512 = vendor.ebx == 0x68747541 and (leaf7.ebx & (1 << 16)) != 0;
        },
        .loongarch64 => {
            const cfg: u64 = loongarch.cpucfg(0x14);
            const mb = ((cfg & 0xffff) + 1) * (@as(u64, 1) << @intCast((cfg >> 16) & 0xff)) * (@as(u64, 1) << @intCast((cfg >> 24) & 0x7f)) / 1024 / 1024;
            c.l3_kb = @intCast(mb * 1024);
            c.prid = loongarch.cpucfg(0);
        },
        else => {},
    }
    return c;
}
