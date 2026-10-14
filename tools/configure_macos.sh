#!/bin/bash
# One-time settings for the macOS test VM (`build/qemutest.zig`), after its
# first login. Serve this directory from the host (QEMU's user network shows
# the host to the guest as 10.0.2.2), e.g. `python -m http.server 8000
# --directory tools`, then in the VM's Terminal:
#   curl -fsSO http://10.0.2.2:8000/configure_macos.sh && bash configure_macos.sh
# Each step says what it does and carries on if it fails; running it again
# is harmless.

step() { echo "== $1"; shift; "$@" || echo "   (failed: $*)"; }

sudo -v || exit 1

# Several settings are guarded by privacy controls (TCC) since Mojave: they
# need Full Disk Access for the app running this (Terminal), whatever `sudo`
# says, and only the user can grant it. The TCC database itself is readable
# only with that access, so it is the test.
if ! sudo cat "/Library/Application Support/com.apple.TCC/TCC.db" > /dev/null 2>&1; then
    echo "Terminal needs Full Disk Access: turn it on in the pane opening now,"
    echo "quit and reopen Terminal, and run this again."
    open "x-apple.systempreferences:com.apple.preference.security?Privacy_AllFiles"
    exit 1
fi

# Verbose boot every time. OpenCore may replace boot-args from its own
# config at boot; if `nvram boot-args` is empty after a reboot, it did.
step "verbose boot" sudo nvram boot-args="-v"

# Boot this volume by default, as General > Startup Disk does; OpenCore
# follows it. With System Integrity Protection on, `bless --setBoot` is not
# permitted (0xe00002e2); then pick the disk in Startup Disk, or select it in
# the OpenCore picker and press Ctrl+Enter.
# if csrutil status | grep -q disabled; then
#     step "default boot disk" sudo bless --mount / --setBoot
# else
#     echo "== default boot disk: skipped, SIP is on (use Startup Disk or Ctrl+Enter in the picker)"
# fi

# Spotlight indexing is the big background CPU user under emulation.
step "spotlight off" sudo mdutil -a -i off

# No sleep, screen saver or disk sleep; a sleeping guest stalls the tests.
step "no sleep" sudo pmset -a sleep 0 displaysleep 0 disksleep 0
step "no screen saver" defaults -currentHost write com.apple.screensaver idleTime -int 0

# Less animation and transparency: the display is drawn in software. Takes
# effect at the next login.
step "reduce motion" defaults write com.apple.universalaccess reduceMotion -bool true
step "reduce transparency" defaults write com.apple.universalaccess reduceTransparency -bool true

# A solid color wallpaper: the dynamic ones flicker without a GPU. Prefers a
# dark gray if this release has one, else the first solid color. macOS asks
# once whether Terminal may control System Events; allow it.
solid="/System/Library/Desktop Pictures/Solid Colors"
picture=$(ls "$solid"/*Gray*.png "$solid"/*.png 2>/dev/null | head -1)
if [ -n "$picture" ]; then
    step "solid wallpaper" osascript -e "tell application \"System Events\" to tell every desktop to set picture to \"$picture\""
else
    echo "== solid wallpaper: no pictures in $solid (use Wallpaper > Colors)"
fi

# Background services the VM does not need. `softwareupdate --schedule` no
# longer changes anything; this preference is what it used to set.
step "siri off" defaults write com.apple.assistant.support "Assistant Enabled" -bool false
step "auto updates off" sudo defaults write /Library/Preferences/com.apple.SoftwareUpdate AutomaticCheckEnabled -bool false
step "time machine off" sudo tmutil disable

# Remote Login (ssh), for `qemu-test-x86_64-macos`.
step "remote login" sudo systemsetup -setremotelogin on

echo "== done; log out and back in for the display settings, and reboot to"
echo "   check verbose boot (nvram boot-args: $(nvram boot-args 2>/dev/null))"
