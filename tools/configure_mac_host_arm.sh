#!/bin/bash
# Settings for an Apple silicon Mac that runs the rigs and the timing over
# ssh (`-Dtest-host`): headless, logged in, never asleep, with macOS's
# background indexing stopped so timings stay steady (src/linalg/blas/
# design.md, Settled decisions: timing on a busy Mac). Run it as the account
# the rigs log in as.
#   bash configure_mac_host_arm.sh          apply what a script can, report the rest
#   bash configure_mac_host_arm.sh --check  only report
# Steps a script cannot take (SIP, iCloud toggles, a password) say what to do
# by hand. Running it again is harmless.

check_only=
[ "$1" = "--check" ] && check_only=1

step() { echo "== $1"; shift; [ -n "$check_only" ] && return 0; "$@" || echo "   (failed: $*)"; }
note() { echo "   $*"; }

# System Integrity Protection must be off first: with it on, launchd drops the
# service overrides below at boot and refuses `bootout` for Apple's services.
# It weakens the machine's protection: only for a headless, physically
# secured machine.
if ! csrutil status 2>/dev/null | grep -q disabled; then
    cat <<'EOF'
System Integrity Protection is on; this script needs it off. Only macOS
Recovery can change it, at the Mac itself:
  1. Shut down. Press and hold the power button until "Loading startup
     options" appears, then choose Options > Continue.
  2. In Recovery, choose Utilities > Terminal and run:  csrutil disable
     (enter the administrator's password when asked), then reboot.
  3. Back in macOS, run this script again.
To turn it back on later, run `csrutil enable` in Recovery the same way.
EOF
    exit 1
fi

[ -z "$check_only" ] && { sudo -v || exit 1; }
uid=$(id -u)

# 1. Logged in at the console. Without a console user `audiomxd` and `configd`
# spin asking for one. Auto-login keeps a session from boot;
# it needs FileVault off and the account password, so it is left to the user:
# System Settings > Users & Groups > Automatically log in as, or over ssh
# `autoLoginUser` in com.apple.loginwindow plus /etc/kcpassword.
echo "== auto-login"
user=$(defaults read /Library/Preferences/com.apple.loginwindow autoLoginUser 2>/dev/null)
if [ -n "$user" ] && [ -f /etc/kcpassword ]; then note "on, as $user"
else note "off: turn it on in Users & Groups (FileVault: $(fdesetup status 2>/dev/null))"; fi

# 2. Never asleep, back after a power cut, reachable on the network. A
# sleeping Mac drops off ssh mid-run.
step "no sleep, wake on network, restart after power loss" \
    sudo pmset -a sleep 0 disksleep 0 standby 0 womp 1 autorestart 1 tcpkeepalive 1
note "$(pmset -g | grep -E '^ (sleep|disksleep|standby|womp|autorestart) ' | tr -s ' ' | tr '\n' ';')"

# 3. Spotlight file indexing off.
step "spotlight indexing off" sudo mdutil -a -i off
note "$(mdutil -s / 2>/dev/null | tail -1 | tr -d '\t')"

# 4. Apple Intelligence off, models removed: RemoveMacAI (github.com/omlahore/
# RemoveMacAI, MIT; its only network call is the app's own Check for Updates).
# Its `off` locks the features through a configuration profile (approved once
# in Device Management) and removes the foundation models; Siri's and speech's
# assets (several GB of UAF sets) stay.
# Missing, it is installed as its README says (the app, from its latest
# release, which `install.sh` checks against the release's SHA-256), then
# turns Apple Intelligence off; the profile it makes still needs approving in
# System Settings > General > Device Management.
echo "== Apple Intelligence (RemoveMacAI)"
# `install.sh app` puts it in /Applications, or ~/Applications if that is not
# writable, and opens it
rma_find() {
    for d in /Applications ~/Applications; do
        [ -x "$d/RemoveMacAI.app/Contents/MacOS/RemoveMacAI" ] && { echo "$d/RemoveMacAI.app/Contents/MacOS/RemoveMacAI"; return; }
    done
}
rma=$(rma_find)
if [ -z "$rma" ] && [ -z "$check_only" ]; then
    note "installing RemoveMacAI"
    curl -fsSL https://raw.githubusercontent.com/omlahore/RemoveMacAI/main/install.sh | bash -s app ||
        note "(failed: install RemoveMacAI by hand, github.com/omlahore/RemoveMacAI)"
    rma=$(rma_find)
    [ -n "$rma" ] && { "$rma" off --yes || note "(failed: $rma off)"; }
    note "approve its profile in System Settings > General > Device Management, then run '$rma off' again"
fi
if [ -n "$rma" ]; then "$rma" status 2>/dev/null | sed -n '1p;/Total/p' | sed 's/^/   /'
else note "RemoveMacAI not installed"; fi

# 5. SIP is off: checked at the top.

# 6. Background services stopped and kept from starting: on-device media,
# text and knowledge indexing, and the Siri AI app (`com.apple.campo`), which
# together can hold several cores long after login. Re-enable one with
# `launchctl enable gui/$UID/<label>` and a reboot.
agents="com.apple.mediaanalysisd com.apple.photoanalysisd com.apple.textunderstandingd
    com.apple.spotlightknowledged com.apple.spotlightknowledged.importer
    com.apple.spotlightknowledged.updater com.apple.corespotlightd
    com.apple.intelligenceplatformd com.apple.campo"
for label in $agents; do
    step "off: $label" sh -c "launchctl bootout gui/$uid/$label 2>/dev/null; launchctl disable gui/$uid/$label"
done
note "disabled: $(launchctl print-disabled gui/$uid | grep -cE "($(echo $agents | tr ' ' '|'))\" => disabled") of $(echo $agents | wc -w)"
# The model managers are system daemons, so `sudo` and the system domain.
daemons="com.apple.modelmanagerd com.apple.modelcatalogd"
for label in $daemons; do
    step "off: $label" sh -c "sudo launchctl bootout system/$label 2>/dev/null; sudo launchctl disable system/$label"
done
note "disabled: $(launchctl print-disabled system | grep -cE "($(echo $daemons | tr ' ' '|'))\" => disabled") of $(echo $daemons | wc -w) system daemons"

# 7. With iCloud signed in, Keychain sync off on this Mac keeps changes to its
# saved Wi-Fi networks local. A setting only: System Settings > your name >
# iCloud > Passwords & Keychain (choose Keep on This Mac).
echo "== iCloud Keychain sync"
db=~/Library/Accounts/Accounts4.sqlite
icloud=$( [ -r "$db" ] && sqlite3 -readonly "$db" "select count(*) from ZACCOUNT a join ZACCOUNTTYPE t
    on a.ZACCOUNTTYPE = t.Z_PK where t.ZACCOUNTTYPEDESCRIPTION = 'iCloud';" 2>/dev/null)
if [ "${icloud:-0}" = 0 ]; then note "no iCloud account signed in"
elif [ -r "$db" ]; then
    on=$(sqlite3 -readonly "$db" "select hex(d.ZNAME) from Z_2ENABLEDDATACLASSES e
        join ZDATACLASS d on e.Z_7ENABLEDDATACLASSES = d.Z_PK
        join ZACCOUNT a on e.Z_2ENABLEDACCOUNTS = a.Z_PK
        join ZACCOUNTTYPE t on a.ZACCOUNTTYPE = t.Z_PK
        where t.ZACCOUNTTYPEDESCRIPTION = 'iCloud';" 2>/dev/null | xxd -r -p | strings | grep -c 'Dataclass.KeychainSync')
    [ "$on" = 0 ] && note "off" || note "on: turn it off in Passwords & Keychain"
else note "no accounts database readable"; fi

# 8. Wi-Fi: the Mac joins its preferred networks in order, so the one it
# should use goes first, and any nearby it must not fall back to can be
# removed: `sudo networksetup -removepreferredwirelessnetwork <device>
# "<name>"` (with Keychain sync off, on this Mac only). Ethernet, or a DHCP
# reservation, keeps the address steady for the ssh alias. Reported only.
echo "== Wi-Fi"
wifi=$(networksetup -listallhardwareports | awk '/Wi-Fi/ {getline; print $2}')
if [ -n "$wifi" ]; then
    nets=$(networksetup -listpreferredwirelessnetworks "$wifi" 2>/dev/null | sed '1d' | tr -d '	')
    note "$(printf '%s
' "$nets" | grep -c .) preferred networks on $wifi; first: $(printf '%s
' "$nets" | head -3 | paste -sd ',' - | sed 's/,/, /g')"
else note "no Wi-Fi device"; fi

# 9. Home, End and Delete in zsh over ssh from other terminals (Windows),
# which send sequences zsh leaves unbound. Added once.
zsh_keys() {
    grep -qs 'Home / End / Delete over SSH' ~/.zshrc && return 0
    cat >> ~/.zshrc <<'EOF'

# Home / End / Delete over SSH from other terminals
bindkey '^[[H'  beginning-of-line   # xterm
bindkey '^[[F'  end-of-line
bindkey '^[OH'  beginning-of-line   # application cursor mode
bindkey '^[OF'  end-of-line
bindkey '^[[1~' beginning-of-line   # VT220 / some Windows terminals
bindkey '^[[4~' end-of-line
bindkey '^[[3~' delete-char         # Delete
EOF
}
step "zsh home and end keys" zsh_keys
grep -qs 'Home / End / Delete over SSH' ~/.zshrc && note "present in ~/.zshrc" || note "missing from ~/.zshrc"

# 10. Remote Login, for the rigs (`-Dtest-host`) and the ssh test steps.
step "remote login" sudo systemsetup -setremotelogin on

echo "== done; reboot, then run with --check to see that the overrides lasted"
