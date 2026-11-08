#!/bin/bash
# Nix as on macmini-m6: the official multi-user install, then nix-darwin
# from this host's directory (its LocalHostName) in nixos-cfg, cloned to
# /etc/nix/git/nixos-cfg and linked as /etc/nix-darwin, darwin-rebuild's
# default, so `sudo darwin-rebuild switch` needs no flags. Run it as the
# account that owns the clone, after configure_mac_host_*.sh, over `ssh -t`:
#   bash configure_mac_nix.sh          set up what is missing
#   bash configure_mac_nix.sh --check  only report
# Needs a key git.syme.biz accepts (~/.ssh/id_ed25519, as macmini-m6's) and
# <host>/flake.nix committed in nixos-cfg. Running it again is harmless.

check_only=
[ "$1" = "--check" ] && check_only=1
host=$(scutil --get LocalHostName)
repo=ssh://git@git.syme.biz:2222/configs/nixos-cfg.git
clone=/etc/nix/git/nixos-cfg
nix=/nix/var/nix/profiles/default/bin/nix
xf=(--extra-experimental-features "nix-command flakes")
# nixpkgs-unstable dropped x86_64-darwin; 26.05 is its last release
case $(uname -m) in
    x86_64) pkgs=github:NixOS/nixpkgs/nixpkgs-26.05-darwin ;;
    *) pkgs=github:NixOS/nixpkgs/nixpkgs-unstable ;;
esac

step() { echo "== $1"; shift; [ -n "$check_only" ] && return 0; "$@" || echo "   (failed: $*)"; }
note() { echo "   $*"; }
stop() { echo "   $*"; exit 1; }

[ -z "$check_only" ] && { sudo -v || exit 1; }

# 1. Nix: Nix Store volume, nixbld users, nix-daemon. NIX_INSTALLER_YES
# answers its prompts; sudo may still ask for the password.
echo "== nix"
if [ -x $nix ]; then note "$($nix --version)"
elif [ -n "$check_only" ]; then stop "not installed"
else
    tmp=$(mktemp -d)
    curl --proto '=https' --tlsv1.2 -sSfL https://nixos.org/nix/install -o "$tmp/install" &&
        NIX_INSTALLER_YES=1 sh "$tmp/install" --daemon || stop "(failed: the Nix installer)"
    rm -rf "$tmp"
fi
. /nix/var/nix/profiles/default/etc/profile.d/nix-daemon.sh

# 2. Nix on PATH for every zsh, non-interactive ssh commands too: nix-darwin's
# /etc/zshenv covers login shells only.
zshenv() {
    grep -qs 'nix-daemon.sh' ~/.zshenv && return 0
    cat >> ~/.zshenv <<'EOF'
# Nix on PATH for every zsh, including non-interactive ssh commands
if [ -e /nix/var/nix/profiles/default/etc/profile.d/nix-daemon.sh ]; then
  . /nix/var/nix/profiles/default/etc/profile.d/nix-daemon.sh
fi
EOF
}
step "nix in ~/.zshenv" zshenv
grep -qs 'nix-daemon.sh' ~/.zshenv && note "present" || note "missing"

# 3. nixos-cfg in /etc/nix/git, owned by this account, cloned with nixpkgs'
# git (no Command Line Tools needed).
echo "== nixos-cfg"
if [ -d $clone/.git ]; then note "cloned at $clone"
elif [ -n "$check_only" ]; then note "not cloned"
else
    [ -f ~/.ssh/id_ed25519 ] || stop "no ~/.ssh/id_ed25519: copy macmini-m6's (git.syme.biz knows it), then run again"
    sudo mkdir -p ${clone%/*} && sudo chown "$USER:staff" ${clone%/*} || stop "(failed: ${clone%/*})"
    GIT_SSH_COMMAND="ssh -o StrictHostKeyChecking=accept-new" \
        $nix "${xf[@]}" shell $pkgs#git -c git clone $repo $clone || stop "(failed: clone $repo)"
fi
[ -n "$check_only" ] && [ ! -d $clone ] && exit 1
# nix sees only files git tracks
[ -f $clone/$host/flake.nix ] ||
    stop "no $host/flake.nix: commit one to nixos-cfg (macmini-m6's is the model), pull it into $clone, run again"

# 4. /etc/nix-darwin, darwin-rebuild's default flake, links to it.
[ -e /etc/nix-darwin ] && [ ! -L /etc/nix-darwin ] && stop "/etc/nix-darwin is not a link: move it aside"
step "/etc/nix-darwin link" sudo ln -sfn $clone/$host /etc/nix-darwin
note "-> $(readlink /etc/nix-darwin)"

# 5. nix-darwin. Built first as this account, so flake.lock is its and a
# broken config stops before /etc changes. Activation refuses /etc files it
# does not know, so they move aside first (macmini-m6 has the same four).
echo "== nix-darwin"
if [ -e /run/current-system/darwin-version.json ]; then
    note "active: $(sed -n 's/.*"darwinLabel": *"\([^"]*\)".*/\1/p' /run/current-system/darwin-version.json); update with 'sudo darwin-rebuild switch'"
elif [ -n "$check_only" ]; then note "not active"
else
    darwin=$(sed -n 's/^ *nix-darwin\.url *= *"\(.*\)";.*/\1/p' $clone/$host/flake.nix)
    [ -n "$darwin" ] || stop "no nix-darwin.url in $host/flake.nix"
    $nix "${xf[@]}" build --no-link "$clone/$host#darwinConfigurations.$host.system" || stop "(failed: build)"
    for f in /etc/nix/nix.conf /etc/bashrc /etc/zshrc /etc/zprofile; do
        if [ -f $f ] && [ ! -L $f ] && [ ! -e $f.before-nix-darwin ]; then
            sudo mv $f $f.before-nix-darwin && note "moved $f to $f.before-nix-darwin"
        fi
    done
    if ! sudo $nix "${xf[@]}" run "$darwin#darwin-rebuild" -- switch --flake "$clone/$host#$host"; then
        note "(failed: switch) fix and run again; to undo the moves: sudo mv <file>.before-nix-darwin <file>"
        exit 1
    fi
    note "active; commit $host/flake.lock from $clone"
fi

echo "== done; open a new shell for nix-darwin's PATH"
