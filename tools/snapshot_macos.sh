qemu-img snapshot -c clean .macos-vm/x86_64-macos/mac_hdd.qcow2
cp .macos-vm/x86_64-macos/OVMF_VARS.fd .macos-vm/x86_64-macos/OVMF_VARS.clean.fd
qemu-img snapshot -l .macos-vm/x86_64-macos/mac_hdd.qcow2
