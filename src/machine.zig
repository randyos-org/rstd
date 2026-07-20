const std = @import("std");
const log = std.log.scoped(.machine);

pub const HardwareInterface = enum {
    acpi,
    dtb,
    none,
};

pub const Firmware = enum {
    bios,
    uefi,
    pftf,
    ofw,
    rpi,
    asahi,
    none,
};

/// firmware-neutral memory region classification; kernel only sees this
/// shape, never raw firmware-native format
pub const MemoryRegionKind = enum {
    usable,
    reserved,
    acpi_reclaimable,
    acpi_nvs,
    bootloader_reclaimable,
    kernel_and_modules,
    mmio,
    bad,
};

/// single physical memory region, firmware-neutral
pub const MemoryRegion = struct {
    phys_start: u64,
    page_count: u64,
    kind: MemoryRegionKind,
};

/// framebuffer pixel channel order, firmware-neutral
pub const FramebufferPixelFormat = enum {
    /// Red, Green, Blue, Reserved
    rgb,
    /// Blue, Green, Red, Reserved
    bgr,
};

/// Video Mode Info
pub const KernelBootVideoModeInfo = struct {
    framebuffer_pointer: [*]volatile u32,
    framebuffer_size: usize,
    horizontal_resolution: u32,
    vertical_resolution: u32,
    pixels_per_scanline: u32,
    pixel_format: FramebufferPixelFormat,
};

/// Kernel Boot Info
pub const KernelBootInfo = struct {
    memory_map: []const MemoryRegion,
    video_mode_info: KernelBootVideoModeInfo,
    hardware_description: ?*anyopaque,
    fw_runtime_ptr: ?*anyopaque,
    kernel_phys_start: usize,
    kernel_phys_end: usize,
    kernel_virt_start: usize,
    kernel_virt_end: usize,
    dwarf_info: *?std.debug.Dwarf,
    /// unix epoch seconds at boot per firmware; null if undetermined
    boot_wall_clock_unix_seconds: ?i64,
};
