#!/bin/bash
# One-time settings for the macOS test VM (`build/qemutest.zig`), after its
# first login. Serve this directory from the host (QEMU's user network shows
# the host to the guest as 10.0.2.2), e.g. `python -m http.server 8000
# --directory tools`, then in the VM's Terminal:
#   curl -fsSO http://10.0.2.2:8000/configure_macos.sh && bash configure_macos.sh
# Each step says what it does and carries on if it fails.

step() { echo "== $1"; shift; "$@" || echo "   (failed: $*)"; }

sudo -v || exit 1

# Verbose boot every time. OpenCore may replace boot-args from its own
# config at boot; if `nvram boot-args` is empty after a reboot, it did.
step "verbose boot" sudo nvram boot-args="-v"

# Boot this volume by default, as General > Startup Disk does; OpenCore
# follows it. Otherwise, select it in the picker and press Ctrl+Enter.
step "default boot disk" sudo bless --mount / --setBoot

# Spotlight indexing is the big background CPU user under emulation.
step "spotlight off" sudo mdutil -a -i off

# No sleep, screen saver or disk sleep; a sleeping guest stalls the tests.
step "no sleep" sudo pmset -a sleep 0 displaysleep 0 disksleep 0
step "no screen saver" defaults -currentHost write com.apple.screensaver idleTime -int 0

# Less animation and transparency: the display is drawn in software. Newer
# releases may ignore these; if so, use Accessibility > Display.
step "reduce motion" defaults write com.apple.universalaccess reduceMotion -bool true
step "reduce transparency" defaults write com.apple.universalaccess reduceTransparency -bool true

# Background services the VM does not need.
step "siri off" defaults write com.apple.assistant.support "Assistant Enabled" -bool false
step "auto updates off" sudo softwareupdate --schedule off
step "time machine off" sudo tmutil disable

# Remote Login (ssh), for `qemu-test-x86_64-macos`. On recent releases this
# may need Full Disk Access for Terminal; if so, use General > Sharing.
step "remote login" sudo systemsetup -setremotelogin on

echo "== done; reboot to check verbose boot (nvram boot-args: $(nvram boot-args 2>/dev/null))"
