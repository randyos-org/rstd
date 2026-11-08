#!/bin/bash
# Settings for an Intel Mac (Sequoia or later) that runs the x86_64 macOS
# tests and timing over ssh: as configure_mac_host_arm.sh, less what Intel
# lacks (Apple Intelligence). Run it as the account the tests log in as.
#   bash configure_mac_host_intel.sh          apply what a script can, report the rest
#   bash configure_mac_host_intel.sh --check  only report
# A fresh Mac first needs Remote Login, which a script cannot turn on without
# Full Disk Access: System Settings > General > Sharing > Remote Login. Then,
# from the PC's MSYS2 shell (alias in ~/.ssh/config):
#   ssh-copy-id macmini-intel
#   scp tools/configure_mac_host_intel.sh macmini-intel:
#   ssh -t macmini-intel bash configure_mac_host_intel.sh
# Steps a script cannot take say what to do by hand. Running it again is
# harmless.

check_only=
[ "$1" = "--check" ] && check_only=1
name=macmini-intel

step() { echo "== $1"; shift; [ -n "$check_only" ] && return 0; "$@" || echo "   (failed: $*)"; }
note() { echo "   $*"; }

[ -z "$check_only" ] && { sudo -v || exit 1; }
uid=$(id -u)

# 1. Remote Login, for the ssh test steps. `systemsetup` needs Full Disk
# Access since Ventura, so it may fail; Sharing always works.
echo "== remote login"
if nc -z 127.0.0.1 22 2>/dev/null; then note "on"
elif [ -z "$check_only" ] && sudo systemsetup -setremotelogin on; then note "turned on"
else note "off: turn it on in System Settings > General > Sharing > Remote Login"; fi

# 2. Host name, for Bonjour (`$name.local`) and the prompt.
set_name() {
    sudo scutil --set ComputerName "$name" && sudo scutil --set LocalHostName "$name" &&
        sudo scutil --set HostName "$name"
}
step "host name $name" set_name
note "$(scutil --get ComputerName) / $(scutil --get LocalHostName).local / $(scutil --get HostName 2>/dev/null)"

# System Integrity Protection must be off for the rest: with it on, launchd
# drops the service overrides below at boot. Only for a headless, physically
# secured machine.
if ! csrutil status 2>/dev/null | grep -q disabled; then
    cat <<'EOF'
System Integrity Protection is on; the rest needs it off. Only macOS
Recovery can change it, at the Mac itself:
  1. Restart, and at once hold Command-R until the Apple logo appears.
  2. In Recovery, choose Utilities > Terminal and run:  csrutil disable
     (enter the administrator's password when asked), then reboot.
  3. Run this script again (over ssh is fine).
To turn it back on later, run `csrutil enable` in Recovery the same way.
EOF
    exit 1
fi

# 3. Logged in at the console, so `audiomxd` and `configd` rest. Auto-login
# needs FileVault off (with it on there is no network before login) and the
# account password: System Settings > Users & Groups > Automatically log in as.
echo "== auto-login"
user=$(defaults read /Library/Preferences/com.apple.loginwindow autoLoginUser 2>/dev/null)
if [ -n "$user" ] && [ -f /etc/kcpassword ]; then note "on, as $user"
else note "off: turn it on in Users & Groups (FileVault: $(fdesetup status 2>/dev/null))"; fi
# Auto-login covers boot only: a wake asks for the password unless the
# screen lock is off. `sysadminctl` prompts for the password.
echo "== screen lock"
lock=$(sysadminctl -screenLock status 2>&1 | sed 's/^.*\] //')
if [ "$lock" = "screenLock is off" ]; then note "off"
elif [ -n "$check_only" ]; then note "on ($lock): a wake asks for the password"
elif sysadminctl -screenLock off -password -; then note "turned off"
else note "(failed: Lock Screen > Require password after screen saver begins or display is turned off > Never)"; fi

# 4. Idle power: display and any hard disk off after 10 minutes (the SSD
# ignores disksleep), sleep after 60. An incoming ssh connection wakes it
# through a Bonjour sleep proxy (a HomePod or Apple TV) or a magic packet; no
# standby or hibernation, so the network stays able to and wakes are quick.
# rstd's ssh tools hold `caffeinate -i` for their runs. Only the settings this
# model has (`pmset -g cap`; no `autopoweroff` on a 2018 Mac mini); `pmset -g`
# can lag, so the saved ones are reported.
power="sleep 60 displaysleep 10 disksleep 10 standby 0 autopoweroff 0 hibernatemode 0 powernap 0 womp 1 autorestart 1 tcpkeepalive 1"
cap=$(pmset -g cap)
args=$(echo $power | xargs -n2 | while read k v; do
    echo "$cap" | grep -qx " $k" && echo "$k $v"; done | xargs)
step "sleep after 60 minutes, wake on network, restart after power loss" sudo pmset -a $args
note "$(pmset -g custom | grep -E '^ (sleep|displaysleep|disksleep|standby|hibernatemode|powernap|womp|autorestart|tcpkeepalive) ' | tr -s ' ' | tr '\n' ';')"

# 5. Spotlight file indexing off.
step "spotlight indexing off" sudo mdutil -a -i off
note "$(mdutil -s / 2>/dev/null | tail -1 | tr -d '\t')"

# 6. Background indexing services stopped and kept from starting; the list
# matches the Apple silicon host's, and a label this release lacks is a
# harmless override. Re-enable one with `launchctl enable gui/$UID/<label>`
# and a reboot. The gui domain exists only while the account is logged in at
# the console (auto-login, step 3). On Sequoia the Spotlight UI agent
# (`com.apple.Spotlight`, Cmd-Space) retries corespotlightd's IndexAgent
# endlessly once that is off (a core in it, half one in launchd), so it goes too.
agents="com.apple.mediaanalysisd com.apple.photoanalysisd com.apple.textunderstandingd
    com.apple.spotlightknowledged com.apple.spotlightknowledged.importer
    com.apple.spotlightknowledged.updater com.apple.corespotlightd
    com.apple.intelligenceplatformd com.apple.campo com.apple.Spotlight"
if launchctl print gui/$uid > /dev/null 2>&1; then
    for label in $agents; do
        step "off: $label" sh -c "launchctl bootout gui/$uid/$label 2>/dev/null; launchctl disable gui/$uid/$label"
    done
    note "disabled: $(launchctl print-disabled gui/$uid | grep -cE "($(echo $agents | tr ' ' '|'))\" => disabled") of $(echo $agents | wc -w)"
else
    echo "== user agents: skipped, no console login; log in at the Mac (or turn on auto-login and reboot), then run again"
fi
daemons="com.apple.modelmanagerd com.apple.modelcatalogd"
for label in $daemons; do
    step "off: $label" sh -c "sudo launchctl bootout system/$label 2>/dev/null; sudo launchctl disable system/$label"
done
note "disabled: $(launchctl print-disabled system | grep -cE "($(echo $daemons | tr ' ' '|'))\" => disabled") of $(echo $daemons | wc -w) system daemons"

# 7. With iCloud signed in, Keychain sync off keeps this Mac's saved Wi-Fi
# networks local: System Settings > your name > iCloud > Passwords & Keychain.
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

# 8. A steady address for the ssh alias: Ethernet, or a DHCP reservation with
# the network's Private Wi-Fi address off (a rotating one breaks the
# reservation): Wi-Fi > the network's Details. Reported only.
echo "== network"
wifi=$(networksetup -listallhardwareports | awk '/Wi-Fi/ {getline; print $2}')
if [ -n "$wifi" ]; then
    nets=$(networksetup -listpreferredwirelessnetworks "$wifi" 2>/dev/null | sed '1d' | tr -d '	')
    note "$(printf '%s\n' "$nets" | grep -c .) preferred networks on $wifi; first: $(printf '%s\n' "$nets" | head -3 | paste -sd ',' - | sed 's/,/, /g')"
else note "no Wi-Fi device"; fi
note "default route via $(route -n get default 2>/dev/null | awk '/interface:/ {print $2}')"

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

echo "== done; reboot, then run with --check to see that the overrides lasted"
