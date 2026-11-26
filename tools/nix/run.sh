#!/usr/bin/env bash
# run.sh <jobs> <qemu argv>...: one rig on the test host (tools/rig/test_host.zig) in its files' directory: waits for a slot (flock in
# /tmp/rstd-qemu-slots), runs QEMU until stdin closes, then removes this and any dead run's directory. `@qemu-data@` is QEMU's data dir.
set -u
d=$PWD
q=
trap 'cd /; rm -rf "$d"' EXIT
trap 'if [ -n "$q" ]; then kill "$q" 2>/dev/null; fi; exit 143' TERM
# a dropped connection must not stop the cleanup
trap '' PIPE
jobs=$1
shift
# `@qemu-data@` in an argument: QEMU's data directory (its firmware)
q_data=$(dirname "$(dirname "$(command -v "$1")")")/share/qemu
q_tag=@qemu-data@
set -- "${@//"$q_tag"/$q_data}"
# this run's lock, held while it lives; a run killed before getting here (in
# the copy) leaves a directory no one holds, which a later run removes a day
# on (older than any run, and than runs from before the lock)
exec 5>"$d/.lock"
flock -n 5
# -H: /tmp is a symlink on macOS
for x in $(find -H /tmp -maxdepth 1 -type d -name 'rstd-qemu.*' -mmin +1440); do
  [ "${x##*/}" != "${d##*/}" ] && flock -n "$x/.lock" true 2>/dev/null && rm -rf "$x"
done
# fd 3: a background job's stdin would be /dev/null
exec 3<&0 </dev/null
{ cat <&3 >/dev/null; kill -TERM $$; } >/dev/null 2>&1 &
w=$!
exec 3<&-

slots=/tmp/rstd-qemu-slots
mkdir -p "$slots"
n=0
while :; do
  for ((i = 0; i < jobs; i++)); do
    exec 4>>"$slots/$i"
    flock -n 4 && break 2
    exec 4>&-
  done
  # about once a minute, which also keeps the build's idle watchdog quiet
  if ((n++ % 30 == 0)); then echo "test-host: waiting for a slot ($jobs busy)"; fi
  sleep 2
done

# QEMU holds fd 4 too, so the slot is free once it has exited
"$@" &
q=$!
wait "$q"
s=$?
kill "$w" 2>/dev/null
exit "$s"
