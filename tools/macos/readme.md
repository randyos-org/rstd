# macOS machines

Setup for the Macs that run the rigs, the tests and the timings over SSH. Each script's header says how to run it.

- `configure_mac_host_arm.sh`: settings for the Apple silicon Mac (`-Dtest-host`, `ssh-test-aarch64-macos`, the timing runs).
- `configure_mac_host_intel.sh`: settings for the Intel Mac (`ssh-test-x86_64-macos`).
- `configure_mac_nix.sh`: Nix on a Mac, after its host script.
