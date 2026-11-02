# Shaderlib

GLSL-friendly math for rgl's shaders (`math.zig`, `sdf.zig`) over `rstd.linalg`
(design.md, "Shaders keep their own API"). This file records what Zig's self-hosted SPIR-V backend
emits for each operation: the native instruction (`OpDot`, the GLSL.std.450 set) or scalar code.
todo.md's S2 item makes the operations native and checks it at build time.

## Probes

`../probe.zig` is one fragment shader per operation (its `Op` enum). `zig build shader-probes
-Dvulkan` builds every probe (`-Dshader-probe=<op>` builds one) and the `rgl`, `glyph` and `image`
shaders, and installs the `.spv` files to `zig-out/spv/`. Disassemble with `spirv-dis` and check
with `spirv-val --target-env vulkan1.2` (Vulkan SDK).

## What the backend emits (S1, 2026-10-30)

0.18.0-dev.751, Debug and ReleaseFast alike (no optimizer runs, so the two differ only in safety
code); every module passes `spirv-val`. Counts are of the probe's function; `s` is scalar, `v3`/`v4`
a vector.

| Operation | Emitted | Native instruction |
| --- | --- | --- |
| `math.dot` (`linalg.vector.dot`) | 3 FMul.s + 2 FAdd.s (scalar `mulAdd` chain) | OpDot |
| `@reduce(.Add, a * b)` | FMul.v3, 3 extracts, 2 FAdd.s | OpDot |
| `dotAsm` (inline asm) | OpDot | OpDot |
| `math.cross` (`vec.cross`) | 6 FMul.s + 3 FSub.s | Cross (68) |
| `crossAsm` (inline asm) | GLSL.std.450 Cross | Cross (68) |
| `length` | scalar dot + Sqrt | Length (66) |
| `distance` | FSub.v3 + scalar dot + Sqrt | Distance (67) |
| `normalize` | scalar dot + Sqrt + FDiv.v3 (+ zero test) | Normalize (69) |
| `reflect` | scalar dot + FMul.v3 + FSub.v3 | Reflect (71) |
| `refract` | scalar dot, Sqrt, vector ops, branch | Refract (72) |
| `faceforward` | scalar dot + FNegate.v3 + branch | FaceForward (70) |
| `inversesqrt` | Sqrt + FDiv (scalar or vector) | InverseSqrt (32) |
| `linalg.vector.mulAdd` | FMul.v3 + FAdd.v3, unfused | Fma (50) |
| `@mulAdd` | GLSL.std.450 Fma | Fma (50) |
| `mix` | 2 FMul + FSub + FAdd, unfused (vector or scalar) | FMix (46) |
| `step_scalar` | branch | Step (48) |
| `smoothstep_scalar` | FMin/FMax + 3 FMul.s, 3 FSub.s, FDiv.s | SmoothStep (49) |
| `fract` | Floor + FSub | Fract (10) |
| `clamp` (`std.math.clamp`) | FMin + FMax + an assert branch | FClamp (43) |
| `@min`/`@max` | FMin/FMax | FMin (37) / FMax (40) |
| `a * b + c` on `@Vector` | FMul.v3 + FAdd.v3 | (vector arithmetic) |
| `matVec` (columns times splats) | 4 FMul.v4 + 3 FAdd.v4 | OpMatrixTimesVector |
| `vecMat` (a dot per column) | 16 FMul.s + 12 FAdd.s | OpVectorTimesMatrix |
| matrix times matrix (`matVec` per column) | FMul.v4/FAdd.v4 per column | OpMatrixTimesMatrix |

- **Native already:** the math builtins (`@sqrt`, `@floor` and the rest), `@min`/`@max` and
  `@mulAdd`. `@Vector` arithmetic stays vector.
- **Never native through rstd or shaderlib:** Dot, Length, Distance, Normalize, Cross, Reflect,
  Refract, FaceForward, InverseSqrt, Fma (rstd's `mulAdd` stays unfused: `fma.zig`'s `fusedMulAddBits`
  is 0 for SPIR-V), FMix, Step, SmoothStep, Fract, FClamp. Each is reachable through inline `asm`
  (`OpDot`; `OpExtInst` with the GLSL.std.450 opcode as an integer), as `dotAsm` and `crossAsm`
  show.
- **Matrices:** shaderlib's `Mat2` to `Mat4` are arrays of columns, and inline asm cannot declare
  `OpTypeMatrix` (issues.md, "SPIR-V inline assembly cannot declare a matrix type"), so the matrix
  instructions are out of reach. Column-major products in vector instructions (`matVec`) are the
  best available; a dot per column (`vecMat`) goes scalar.
- **Branches:** every inlined call keeps its structured blocks (phi-carried returns), even in
  ReleaseFast (44 in `vecMat`); `spirv-opt` would fold them.
- **Fixed on the way:** shaderlib's `inversesqrt`, `reflect` and `refract` did not compile for
  vectors (`2 * dot(i, n) * n` mixes a scalar and a vector); they now splat their scalars.
- **Compiler crash:** `std.fmt.comptimePrint` as an `asm` template crashes the backend (issues.md),
  so templates are literal.
