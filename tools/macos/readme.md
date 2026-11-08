# macOS machines

The macOS VM rig (`build/qemu/macos.zig`) and the Apple silicon Mac that runs the rigs and timings over SSH.

- `configure_macos_vm.sh`: one-time settings inside the VM (its header says how to run it).
- `configure_mac_host_arm.sh`: settings for the Mac host (`-Dtest-host`, the timing runs).
- `install_macos.sh`, `snapshot_macos.sh`: by-hand helpers for the VM's install media and a snapshot.

## The macOS VM

rstd's unit tests in an installed macOS, which the Mach-O kernels need, copied in and run over SSH. It boots with OSX-KVM's OpenCore and firmware, a lazy dependency fetched only when rstd is the project being built (`-Dmacos-vm-off` skips it). The VM lives in `.macos-vm/<name>/` (gitignored).

Set it up once:

1. `zig build qemu-setup-<name>`: fetches a recovery image (`-Dmac-release`, Sonoma by default) and makes the disk.
2. `zig build qemu-install-<name>`, then in the window: "macOS Base System" in the OpenCore picker; Disk Utility, erase the QEMU HARDDISK as APFS; reinstall macOS onto it (pick "macOS Installer", then the disk, at each reboot); create an account; run `tools/macos/configure_macos_vm.sh` in it, which makes the disk the default boot, turns on Remote Login and turns off what slows an emulated Mac.
3. Authorize an SSH key (`ssh-keygen` makes one if you have none): `type %USERPROFILE%\.ssh\id_ed25519.pub | ssh -p 2222 <account>@localhost "mkdir -p ~/.ssh && cat >> ~/.ssh/authorized_keys"`.

Then `zig build qemu-start-<name>` boots it in a window, and, while it is up, `zig build qemu-test-<name> -Dmac-user=<account>` runs the tests in it. With the VM shut down, `qemu-snapshot-<name>` saves a restore point (both disks and the NVRAM) and `qemu-restore-<name>` goes back to it (`-Dmac-snapshot=<name>`, default `clean`): worth doing once it is set up, since a hard stop of the VM can leave macOS unable to boot.
