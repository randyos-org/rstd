#!/usr/bin/env bash
# One qemu rig on the test host (tools/test_host.zig), run in the directory
# holding its files, inside the flake's shell. Waits for one of <jobs> slots
# (flock(1) locks in /tmp/rstd-qemu-slots, so across builds too; a dead
# holder's lock goes with it), runs QEMU, ends it when stdin closes (the
# build's ssh is gone), and removes the directory.
#
#   run.sh <jobs> <qemu argv>...
set -u
d=$PWD
q=
trap 'cd /; rm -rf "$d"' EXIT
trap 'if [ -n "$q" ]; then kill "$q" 2>/dev/null; fi; exit 143' TERM
# a dropped connection must not stop the cleanup
trap '' PIPE
jobs=$1
shift
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
