# Notes for agents working in rstd

Rules and lessons for agents. The RandyOS repo's `CLAUDE.md` applies too when rstd is checked out inside it. Notes about one machine (paths, hosts, caches) go in `CLAUDE.local.md`, which is not committed. Brief background agents to read both first.

## The backlog: GitHub issues

- Open work is the [issues](https://github.com/randyos-org/rstd/issues) of `randyos-org/rstd`. Each title starts with its ID (`N37: ...`); the docs cite IDs, and a search for the ID finds the issue. Within issues, refer to another as `#N`.
- The current milestone is the task list in the open `dev` to `main` pull request's description, in work order, top to bottom. Check an issue off there when it lands on `dev`, and close the issue with a comment naming the commit.
- Every issue is one atomic task, landed in one commit: no multistep issues. Split a plan into one issue per step, each with its own ID (user, 2026-11-02).
- IDs never change. A new issue takes the next number of its prefix: M, S, T and Z are linear algebra's, H housekeeping's, N the rest's. At the move from `todo.md` (2026-11-08) the highest were H149, M18, N52, S2, T164 and Z6. A split issue's ID is retired: M3 became M8 to M16, and older docs and commits still say M3.
- `**Needs:**` links the issues that must land first; `**Blocked:**` names an outside cause (hardware, an upstream fix, a decision), with the `blocked` label.
- Labels: `backlog` for what is outside the milestone (tier 3 and 4 bugs, follow-ups, tuning, optional and blocked work); `awaiting input` for a decision the user makes, mentioning `@emanspeaks`; `housekeeping` for readability and structure work; `blocked`. Pull requests take the release-note labels of `.github/release.yml` (`fix`, `feature`, `performance`, `tests`, `refactor`, `ci`, `skip-changelog`).
- Tiers ([`support.md`](support.md)): tier 1 and 2 bugs go in the milestone; the rest are `backlog`. Don't open issues for other architectures, follow-ups or tuning on your own; report them.
- Background agents don't edit issues or the pull request's list; they report what to check off or open.

## Git

- One agent, one task at a time, with the user. Each task gets a branch from `dev` (`git switch -c wip/<ID> dev`, after checking `git status` for the user's edits); verified pieces commit there. It lands on `dev` only as the user says, then the branch is deleted. Unfinished work stays on its branch or as unstaged edits, and its issue says what is left, so nothing lives only in a chat or a scratchpad.
- `main` takes only releases, through the `dev` to `main` pull request; its ruleset requires a pull request and passing checks. CI ([`.github/workflows/ci.yml`](.github/workflows/ci.yml)) checks pull requests and `main`, tags each release `v<version>` from `build.zig.zon`, and bumps `dev` to the next patch version. A pull request to `main` that changes the package needs a version newer than the last tag.
- Never stage the user's edits, push, or touch `main` unless asked: the user reviews diffs.
- The user stages their own edits in the trunk. Before any commit, check `git diff --cached` and `git diff --name-only` against your paths, so a commit takes nothing of theirs; build each commit from exactly your hunks.
- Commit messages are one short line (about 70 characters), with no message body. No `Co-Authored-By` or any AI attribution, on any branch. Commits are the user's (Randy Eckman), `--no-gpg-sign`.
- `CLAUDE.md` is committed; `CLAUDE.local.md` never is.
- Several agents may share the trunk: touch only your own hunks, and say so before changing another agent's.
- After a rejected or interrupted command, check what it already changed before going on; don't regenerate or edit `src/` while the user has a build running.
- GitHub is canonical; [`.github/workflows/mirror-codeberg.yml`](.github/workflows/mirror-codeberg.yml) mirrors every branch and tag to Codeberg for dependents still fetching from there. Never push to Codeberg: the mirror overwrites it.

## Design rules

- `rstd.linalg` is an anthology (user, 2026-10-30): for each routine, the best open-source method (OpenBLAS, BLIS, ReLAPACK, faer, Reference-LAPACK and so on), ported as an algorithm into Zig, not line for line. A Zig-native, column-major BLAS: no CBLAS or Fortran names, no `xerbla`, no row-major layer.
- Port faithfully: each ported function's doc comment names its source file; thresholds and parameters come from the reference, cited, never magic numbers. Write complete code: implementation, tests, docs and generator input; no stubs left behind unless the issue allows it.
- Reproducibility beats speed: for a given core, results are bit-identical for any worker count, `Io` and scheduling. Reductions combine in a fixed, size-determined order. A vendor backend (Accelerate) is verified reproducible or runs single-threaded.
- OpenBLAS is the battle-tested baseline: replacing or dropping one of its kernels needs a real defect, benchmarks on the affected cores and the user's OK. Patch the converted assembly first.
- Licenses: MIT is fine. BSD-3 (OpenBLAS, BLIS, PocketFFT) and pffft's UCAR license go in with a `LICENSE` section in the same commit as the first ported code, plus a provenance line in each ported file; a ported header naming a holder `LICENSE` lacks adds it in that commit. Anything else: stop and ask.
- Kernels are reached by slot only (user, 2026-11-07): `portable.fetch(T, op, getter)` in generic code, `blas.kernel(.slot)` for a fixed one; no named per-slot wrappers. New code doesn't build slot names from strings (H115 to H119 convert the existing ones).
- New APIs: panics only for unrecoverable states (a broken invariant, possible memory corruption from assembly); bad input is an error. Exact error sets with plain names (`SliceTooShort`, not `InvalidLda`); slices, not pointers plus sizes. The existing BLAS gets this in Z1 to Z5, not piecemeal.
- Target platforms come first: RandyOS's (x86_64, aarch64, 32-bit arm, PowerPC 750, UEFI) and zfls's (Windows, Linux, macOS on x86_64 and aarch64) must work and be tested. Elsewhere, where nothing can test it but OpenBLAS and BLIS agree, follow them and document it; record failures rather than sinking hours into them.

## Zig style

- `const sysinfo = @import("builtin");`, never `builtin`. Every `@import` is a `const` near the top of the file (except per-arch switches, generated role tables and test code).
- A switch mapping tags to values has one prong per distinct value (`.power6, .power8 => 128,`). Fix generated code in its generator (`tools/blas_asm`), never by hand.
- Files ideally under 500 lines, one per logical unit (a table, a routine, a concern), a folder for a unit too complex for one file; never merge distinct units to cut the file count, only fold thin glue into its home; code moved out of a big file gets its own file (user, 2026-11-07). One module list and `inline for` over it, never a handwritten chain; string identifiers parse into an enum and an exhaustive `switch`; code shared by routines goes in `util.zig`, never copied.
- Allocator first; `try`, never `catch @panic("OOM")`; a fn-pointer type's error set names the real set (`std.mem.Allocator.Error`).
- `src/_os/` keeps its underscore: platform-specific code that the public wrappers (`io.zig` and the like) point to, not an API of its own (user, 2026-11-07).
- Test-only files go in `_test/` folders; moving one renames its tests, so update the rig filters in `build/qemu/machines.zig`.
- Keep `buildutils.OptimizeMode`. Comments are sparse and terse, no em dashes (semicolons instead), American English. Don't remove comments or commented-out lines you didn't write.
- Data derived from an outside source says where from (file, commit) and that it is regenerated, not hand-edited. No Python in the repo: helper scripts stay in scratch and are rewritten when needed.
- Markdown: one line per paragraph or bullet (the editor wraps); lint clean with `npx --no-install markdownlint-cli2 <files>`. Don't add dictionary words; list new ones in your report.

## Toolchain traps

- The user's local Zig is an experimental build (0.18 dev with the fork-only `@Matrix`). A compiler crash or odd codegen there may be its bug: stop and report it with the traceback; never work around it in rstd.
- `@Matrix` exists only in that build: it lives in `src/linalg/_matrix/native.zig`, chosen as the `rstd_matrix` module when the compiler has it, else `fallback.zig`. Zfls builds rstd on Zig 0.17, as CI does: build-script changes must pass `zig build -h` on 0.18 and 0.17.
- `**` array repetition is gone; use `@splat`. The error only says "whitespace on one side".
- `zig fmt` may rewrite `@intFromEnum` to `@backingInt` on backed enums: reread the diff and retest after formatting.
- Zig 0.18 exits 0 on inline-asm errors, silently dropping the line: read the output for `error: <inline asm>` even at exit 0.
- A build killed while linking leaves a cache entry with no executable, and later builds fail with `FileNotFound`: report it, or move that entry's manifest (`h/`) and output (`o/`) aside; never wipe or edit a shared cache.
- Never edit, move or delete anything in a `zig-pkg/`; fix a dependency by repinning `build.zig.zon` to a fixed upstream or the user's fork.
- Code behind an OS or arch switch compiles only for that target: compile-check every target it serves (at least `aarch64-linux-gnu`, `x86_64-linux-gnu`, `aarch64-windows`, `aarch64-macos`). Code every rig compiles must build for 32-bit targets without 64-bit atomics (`qemu-check-arm-softfp`, `qemu-check-powerpc-750`).

## Testing and runs

- Cheap checks per change: native filtered tests (`zig build rstd-test -Drstd-test-filter=<name>`, repeatable; `-Dblas-cores=all` for every core the host runs), `qemu-check-<rig>` for touched arches (never the all-rigs `qemu-check`, up to 80 minutes), `zig build -h` on 0.18 and 0.17 for build changes. Test only what a change reaches. Rig and Mac runs come in batches after a section of work.
- Foreign targets run on the QEMU rigs (`build/qemu/`; `zig build -h` lists them), never WSL, and never in the trunk. Never run a state-changing step (the macOS VM's snapshot or install steps) to check that it compiles; `zig build -h` does that.
- A rig's live log is `.zig-cache/qemu-test-<rig>.log` in its tree. Rigs print `[core]` as each core starts; a watchdog ends a rig after 60 minutes without output. In debug builds every `inline for` iteration keeps its own frame, so per-core loops call a helper per iteration.
- A backgrounded `(build; echo exit $?)` always reports success: read the log's own exit line and `: error:` lines.
- Interrupt runs of code about to change. Kill only processes you started, matched by their own command line; never `taskkill /T` or a tree kill from a parent PID (Windows reuses PIDs), and never an `ssh` session or anything interactive.
- Never write a redirect that can create a `nul` file: `/dev/null` in Bash, `$null` in PowerShell.
- Search with `rg` (with an explicit path, or it reads stdin and hangs) or the Grep tool, excluding `zig-pkg`, `vendor`, `.zig-cache` and `.agent`; recursive `grep -r` crawls the caches for minutes.
- No ad hoc tool installs (markdownlint excepted); use what is on PATH. Prefer WebFetch or curl; a browser opens windows on the user's desktop, so close every page you opened and say so.

## Upstream problems

- Each problem found in another project (Zig, LLVM, QEMU, OpenBLAS and so on) is an issue here as it is found, labeled `upstream: <project>` (`upstream: zig`, `upstream: llvm`, ...) and with ISA labels where it is specific: what we saw, the workaround and where, the text of an upstream report. Never file it upstream; the user triages them and links the upstream report in the issue once filed.
- Code and docs cite the issue, `issue #186` in comments and a link in Markdown, not a file of notes; drop the workaround and close the issue once upstream fixes it.
- T52 ([#139](https://github.com/randyos-org/rstd/issues/139)) lists the workarounds to recheck on each Zig update.
