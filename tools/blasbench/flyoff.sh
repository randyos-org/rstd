#!/usr/bin/env bash
# `flyoff.sh zen3 <OpenBLAS dll>` or `flyoff.sh m6 [--collect]`: the flyoff's timing runs under the timing rule (readme.md), from the
# rstd root under MSYS2's bash ($ZIG, $MAC default macmini-m6). `flyoff.sh run <macos|windows> <dir> <name>:<sections>:<command>...`
# is each host's engine: <command> in <dir> with `@` as the sections, output to runs/<name>-<i>.log. From the R4 runs' script (2026-10-29).
set -u
here=$(cd "$(dirname "$0")" && pwd)
zig=${ZIG:-zig}
mac=${MAC:-macmini-m6}
say() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a "$D/progress.log"; }
TAB=$(printf '\t')

# ---- the idle rule, per OS

# one idle check, a line in poll.log; 0 when idle
idle_macos() {
  local q l out procs top5 maxp lim why
  q=$(ps -Ao comm= | grep -c qemu-system)
  l=$(sysctl -n vm.loadavg | awk '{print $2}')
  out=$(top -l 2 -s 5 -o cpu -n 14 -stats pid,command,cpu | awk '/^PID/{n++; next} n==2')
  procs=$(echo "$out" | awk 'NF>=3 {cpu=$NF; $1=""; $NF=""; sub(/^ +/,""); sub(/ +$/,""); print cpu "\t" $0}' | awk -F '\t' '$2!="top" && $2 !~ /^sshd/' | sort -t "$TAB" -k1,1 -g -r)
  top5=$(echo "$procs" | head -5 | awk -F '\t' '{printf "%s %s; ", $2, $1}')
  # the user keeps Activity Monitor open: it and sysmond are out of the 5%
  # test, their share added to the load limit
  maxp=$(echo "$procs" | awk -F '\t' '$2 !~ /^(Activity Monitor|sysmond)$/' | head -1 | cut -f1)
  lim=$(echo "$procs" | awk -F '\t' '$2 ~ /^(Activity Monitor|sysmond)$/ {s+=$1} END{printf "%.2f", 2+s/100}')
  why=""
  [ "$q" != 0 ] && why="${why}qemu "
  awk -v l="$l" -v m="$lim" 'BEGIN{exit !(l>=m)}' && why="${why}load "
  awk -v t="${maxp:-0}" 'BEGIN{exit !(t>5)}' && why="${why}proc "
  local st=IDLE
  [ -n "$why" ] && st="BUSY(${why% })"
  echo "$(date '+%Y-%m-%d %H:%M:%S') load=$l (under $lim) qemu=$q top5=[${top5%; }] $st $1" >> "$D/poll.log"
  [ -z "$why" ]
}

# Windows' per-process CPU (% of one core) and processor queue over 5 s, as
# "<cpu>\t<name>" lines and a last "queue\t<length>" line
win_sample() {
  typeperf '\Process(*)\% Processor Time' '\System\Processor Queue Length' -si 5 -sc 2 2>/dev/null | tr -d '\r' | awk -F '","' '
    /^"\(PDH-CSV/ { for (i = 2; i <= NF; i++) { h = $i; gsub(/"/, "", h); if (h ~ /Processor Queue/) name[i] = "queue"; else { sub(/.*Process\(/, "", h); sub(/\).*/, "", h); name[i] = h } } next }
    /^"[0-9]/ { last = $0 }
    END { n = split(last, v, "\",\""); for (i = 2; i <= n; i++) { x = v[i]; gsub(/"/, "", x); if (name[i] == "queue") q = x; else if (name[i] != "_Total" && name[i] != "Idle") printf "%.1f\t%s\n", x, name[i] } printf "queue\t%s\n", q }'
}

# the same rule here: the processor queue under 2 stands for the load; this
# script's own tools are out of the 5% test
idle_windows() {
  local q s procs top5 maxp ql why
  q=$(tasklist //FO CSV //NH 2>/dev/null | grep -ci '^"qemu-system')
  s=$(win_sample)
  ql=$(echo "$s" | awk -F '\t' '$1=="queue"{print $2}')
  procs=$(echo "$s" | awk -F '\t' '$1!="queue" && $2 !~ /^(typeperf|bash|sh|ssh|conhost|tasklist|awk|Taskmgr)(#[0-9]+)?$/' | sort -t "$TAB" -k1,1 -g -r)
  top5=$(echo "$procs" | head -5 | awk -F '\t' '{printf "%s %s; ", $2, $1}')
  maxp=$(echo "$procs" | head -1 | cut -f1)
  why=""
  [ "$q" != 0 ] && why="${why}qemu "
  awk -v l="${ql:-0}" 'BEGIN{exit !(l>=2)}' && why="${why}queue "
  awk -v t="${maxp:-0}" 'BEGIN{exit !(t>5)}' && why="${why}proc "
  local st=IDLE
  [ -n "$why" ] && st="BUSY(${why% })"
  echo "$(date '+%Y-%m-%d %H:%M:%S') queue=$ql qemu=$q top5=[${top5%; }] $st $1" >> "$D/poll.log"
  [ -z "$why" ]
}

waitidle() { until "idle_$rule" "$1"; do sleep 30; done; }

snap() {
  echo "$1 $(date '+%Y-%m-%d %H:%M:%S %Z')"
  if [ "$rule" = macos ]; then
    uptime
    top -l 2 -s 5 -n 8 -o cpu -stats command,cpu | tail -9
    echo "QEMU $(ps -Ao comm= | grep -c qemu-system)"
  else
    win_sample | sort -t "$TAB" -k1,1 -g -r | head -9
    echo "QEMU $(tasklist //FO CSV //NH 2>/dev/null | grep -ci '^"qemu-system')"
  fi
}

# one line in the watch file when another process is over 5% or another
# QEMU runs; $1 the QEMUs the run itself starts
watch1() {
  if [ "$rule" = macos ]; then
    ps -Ao pcpu=,comm= | awk -v own="$1" -v t="$(date +%H:%M:%S)" '{c=$2; sub(/.*\//,"",c)} c ~ /^qemu-system/ {q++} $1>5 && c !~ /^(blas_bench|qemu-system|sshd|nix|Activity|sysmond|ps$)/ {s=s c "=" $1 " "} END {if (q>own) s=s "qemu=" q; if (s!="") print t, s}'
  else
    win_sample | awk -F '\t' -v own="$1" -v t="$(date +%H:%M:%S)" '$1=="queue"{next} $2 ~ /^qemu-system/ {q++} $1>5 && $2 !~ /^(blas_bench|qemu-system|typeperf|bash|sh|ssh|conhost|tasklist|awk|Taskmgr)(#[0-9]+)?$/ {s=s $2 "=" $1 " "} END {if (q>own) s=s "qemu=" q; if (s!="") print t, s}'
  fi
}

# one run: name sections command
one() {
  local name=$1 secs=$2 cmd=$3 t pid wpid kpid rc start own=0
  case $cmd in *qemu-system*) own=1 ;; esac
  for t in 1 2 3 4; do
    waitidle "(before $name)"
    say "start $name try $t: $secs"
    snap BEFORE > "$R/$name.before" 2>&1
    : > "$R/$name.watch"
    (while :; do
      watch1 "$own" >> "$R/$name.watch"
      sleep 20
    done) &
    wpid=$!
    start=$(date +%s)
    (cd "$D" && eval "${cmd//@/$secs}") > "$R/$name.log" 2>&1 < /dev/null &
    pid=$!
    (sleep 7200 && kill $pid) 2> /dev/null &
    kpid=$!
    wait $pid
    rc=$?
    kill $kpid $wpid 2> /dev/null
    wait $kpid $wpid 2> /dev/null
    {
      snap AFTER
      echo "rc=$rc secs=$(($(date +%s) - start))"
    } > "$R/$name.after" 2>&1
    say "end $name rc=$rc watch=$(wc -l < "$R/$name.watch" | tr -d ' ')"
    if [ ! -s "$R/$name.watch" ] || [ $t = 4 ]; then return 0; fi
    for f in before after log watch; do mv "$R/$name.$f" "$R/$name.dist$t.$f"; done
    say "disturbed $name: $(tr '\n' ' ' < "$R/$name.dist$t.watch")"
  done
}

# a configuration twice, then its noisy sections up to 3 more times
cfg() {
  local c=$1 all=$2 cmd=$3 i secs
  for i in 1 2; do one "$c-$i" "$all" "$cmd"; done
  for i in 3 4 5; do
    "$merge" --merge "$R/$c.merged.md" $(ls "$R/$c"-[0-9].log | sort) > "$R/$c.report"
    secs=$(grep '^NOISY' "$R/$c.report" | awk '{print $2}' | sort -u | tr '\n' ' ')
    say "$c: $(grep -c '^NOISY' "$R/$c.report") noisy rows before run $i [$secs]"
    [ -z "$secs" ] && break
    one "$c-$i" "$secs" "$cmd"
  done
  "$merge" --merge "$R/$c.merged.md" $(ls "$R/$c"-[0-9].log | sort) > "$R/$c.report"
  say "$c final: $(grep -c '^NOISY' "$R/$c.report") noisy, $(grep -c MISMATCH "$R/$c.report") mismatches"
}

# the engine: rule dir config...
engine() {
  rule=$1
  D=$(cd "$2" && pwd)
  shift 2
  R=$D/runs
  mkdir -p "$R"
  merge=$D/blas_bench
  [ -x "$merge.exe" ] && merge=$merge.exe
  trap 'say exit; [ -e "$D/DONE" ] || echo "ended early $(date "+%F %T"); see progress.log" > "$D/DONE"' EXIT
  rm -f "$D/DONE"
  # the first window: idle, then idle again 2 minutes on
  while :; do
    waitidle "(first window)"
    sleep 120
    "idle_$rule" "(first window, recheck)" && break
  done
  say "idle window"
  local spec name secs cmd
  for spec in "$@"; do
    name=${spec%%:*}
    secs=${spec#*:}
    cmd=${secs#*:}
    secs=${secs%%:*}
    cfg "$name" "$secs" "$cmd"
  done
  say ALLDONE
  echo "finished $(date '+%F %T'); see progress.log and runs/" > "$D/DONE"
}

# the VM's command for the engine: QEMU's arguments (vm.args), the kernel and
# initramfs, the bench's arguments (vm.append), threads $2 and `@`; $1 runs
# QEMU's argv (Nix's shell on the Mac)
vm_cmd() {
  local args
  args=$(tr '\n' ' ' < "$D/vm.args")
  echo "$1 $args -kernel linux -initrd initramfs.cpio -append \"$(cat "$D/vm.append") --threads=$2 @\" | tr -d '\\r'"
}

# the one-thread sections, and the threaded ones (indicative)
all_secs="lib1 lib2 lib3 lapack"
thr_secs="lib3 lapack"

case ${1:-} in
run)
  shift
  engine "$@"
  ;;
zen3)
  [ $# -ge 2 ] || {
    echo "usage: flyoff.sh zen3 <OpenBLAS dll>" >&2
    exit 2
  }
  dll=$(cygpath -w "$2" 2> /dev/null || echo "$2")
  "$zig" build flyoff-stage-zen3 "-Dtest-host=$mac" || exit 1
  D=$(pwd)/zig-out/flyoff-zen3
  bench="./blas_bench.exe --openblas='$dll'"
  engine windows "$D" \
    "win-1t:$all_secs:$bench --threads=1 @" \
    "win-thr:$thr_secs:$bench --threads=default @" \
    "lx-1t:$all_secs:$(vm_cmd '' 1)" \
    "lx-thr:$thr_secs:$(vm_cmd '' default)"
  ;;
mac-run)
  # on the Mac, in the staged directory: the libraries, then the runs
  D=$(cd "$2" && pwd)
  export PATH=/nix/var/nix/profiles/default/bin:$PATH
  rev=$(sed -n 's/.*"rev": *"\([0-9a-f]*\)".*/\1/p' "$D/nix/flake.lock" | head -1)
  : > "$D/progress.log"
  say "flyoff start; macOS $(sw_vers -productVersion) ($(sw_vers -buildVersion))"
  bash "$D/nix/blas.sh" aarch64-darwin "$rev" "$D/bench-nix" >> "$D/progress.log" 2>&1 || {
    say "blas.sh darwin failed"
    echo "ended early: blas.sh darwin" > "$D/DONE"
    exit 1
  }
  bash "$D/nix/blas.sh" aarch64-linux "$rev" "$D/bench-nix" > "$D/libs.cpio" 2>> "$D/progress.log" || {
    say "blas.sh linux failed"
    echo "ended early: blas.sh linux" > "$D/DONE"
    exit 1
  }
  cat "$D/base.cpio" "$D/libs.cpio" > "$D/initramfs.cpio"
  L=/tmp/rstd-blas/aarch64-darwin
  bench="./blas_bench --openblas=$L/openblasCompat/lib/libopenblas.dylib --reference=$L/lapack-reference/lib/liblapack.dylib --eigen=$L/eigen/lib/libeigen_shim.dylib"
  qemu="nix shell ./nix#qemu -c"
  engine macos "$D" \
    "mac-accel:$all_secs:$bench --threads=1 @" \
    "mac-kern:$all_secs:$bench --threads=1 --rstd-accelerate=0 @" \
    "mac-thr:$thr_secs:$bench --threads=default @" \
    "lx:$all_secs:$(vm_cmd "$qemu" 1)" \
    "lx-thr:$thr_secs:$(vm_cmd "$qemu" default)"
  ;;
m6)
  D=$(pwd)/zig-out/flyoff-m6
  ssh=(ssh -o BatchMode=yes -o ConnectTimeout=30 -o ServerAliveInterval=30 "$mac")
  if [ "${2:-}" != --collect ]; then
    "$zig" build flyoff-stage-m6 || exit 1
    rm -rf "$D/runs"
    "${ssh[@]}" 'rm -rf rstd-flyoff' || exit 1
    scp -q -r "$D" "$mac:rstd-flyoff" || exit 1
    "${ssh[@]}" 'cd rstd-flyoff && chmod +x blas_bench && (nohup caffeinate -i bash flyoff.sh mac-run . > run.log 2>&1 < /dev/null &)' || exit 1
    echo "flyoff: started on $mac in ~/rstd-flyoff; waiting for DONE (Ctrl-C leaves it running; collect with flyoff.sh m6 --collect)"
    until "${ssh[@]}" 'test -e rstd-flyoff/DONE' 2> /dev/null; do sleep 300; done
  fi
  "${ssh[@]}" 'cat rstd-flyoff/DONE'
  scp -q -r "$mac:rstd-flyoff/runs" "$mac:rstd-flyoff/progress.log" "$mac:rstd-flyoff/poll.log" "$D/" && echo "flyoff: results in $D/runs"
  ;;
*)
  sed -n '2,27p' "$0" | sed 's/^# \{0,1\}//'
  exit 2
  ;;
esac
