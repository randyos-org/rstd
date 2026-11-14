# Shaderlib

GLSL-friendly math for rgl's shaders (`math.zig`, `sdf.zig`) over `rstd.linalg` (design.md, "Shaders keep their own API"). This file records what Zig's self-hosted SPIR-V backend emits for each operation: the native instruction (`OpDot`, the GLSL.std.450 set) or scalar code.

## Probes

`../probe.zig` is one fragment shader per operation (its `Op` enum). `zig build shader-probes -Dvulkan` builds every probe (`-Dshader-probe=<op>` builds one) and the `rgl`, `glyph` and `image` shaders, and installs the `.spv` files to `zig-out/spv/`. Each probe with a native instruction (`probe.expect`) is disassembled with `spirv-dis` (Vulkan SDK, on `PATH`), and the step fails if the instruction is missing, so an operation that falls back to scalar code breaks the build. Check modules with `spirv-val --target-env vulkan1.2`.

## What the backend emits (S1, 2026-10-30)

0.18.0-dev.751, Debug and ReleaseFast alike (no optimizer runs, so the two differ only in safety code); every module passes `spirv-val`. Counts are of the probe's function; `s` is scalar, `v3`/`v4` a vector.

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
| `stepScalar` | branch | Step (48) |
| `smoothstepScalar` | FMin/FMax + 3 FMul.s, 3 FSub.s, FDiv.s | SmoothStep (49) |
| `fract` | Floor + FSub | Fract (10) |
| `clamp` (`std.math.clamp`) | FMin + FMax + an assert branch | FClamp (43) |
| `@min`/`@max` | FMin/FMax | FMin (37) / FMax (40) |
| `a * b + c` on `@Vector` | FMul.v3 + FAdd.v3 | (vector arithmetic) |
| `matrix.mul(m, v)` | OpMatrixTimesVector | OpMatrixTimesVector |
| `matrix.mul(v, m)` | OpVectorTimesMatrix | OpVectorTimesMatrix |
| `matrix.mul(m, n)` | OpMatrixTimesMatrix | OpMatrixTimesMatrix |

- **Native already:** the math builtins (`@sqrt`, `@floor` and the rest), `@min`/`@max` and `@mulAdd`. `@Vector` arithmetic stays vector.
- **Never native through rstd or shaderlib (at S1):** Dot, Length, Distance, Normalize, Cross, Reflect, Refract, FaceForward, InverseSqrt, Fma (rstd's `mulAdd` stays unfused: `fma.zig`'s `fusedMulAddBits` is 0 for SPIR-V), FMix, Step, SmoothStep, Fract, FClamp. Each is reachable through inline `asm` (`OpDot`; `OpExtInst` with the GLSL.std.450 opcode as an integer), as `dotAsm` and `crossAsm` show. M12 made `linalg`'s native (below).
- **Matrices**: shaderlib's `Mat2` to `Mat4` are `linalg.matrix.StaticMat(N, N, f32)`, whose products go through the fork's `@Matrix`, so they are SPIR-V's matrix instructions. Loading the packed `[C][R]f32` into `OpTypeMatrix` costs `OpCompositeExtract`s (36 for a 4 by 4). Before, arrays of columns gave 4 FMul.v4 and 3 FAdd.v4 for `m v` and 16 FMul.s and 12 FAdd.s for `v m`; inline asm still cannot declare `OpTypeMatrix` ([#198](https://github.com/randyos-org/rstd/issues/198)).
- **Branches:** every inlined call keeps its structured blocks (phi-carried returns), even in ReleaseFast (44 in `vecMat`); `spirv-opt` would fold them.
- **Fixed on the way:** shaderlib's `inversesqrt`, `reflect` and `refract` did not compile for vectors (`2 * dot(i, n) * n` mixes a scalar and a vector); they now splat their scalars.
- **Compiler crash:** `std.fmt.comptimePrint` as an `asm` template crashes the backend ([#199](https://github.com/randyos-org/rstd/issues/199)), so templates are literal.

## `linalg`'s SPIR-V forms (M12, 2026-11-07)

0.18.0-dev.859; every module passes `spirv-val`. `linalg.vector` on a `@Vector` of 2 to 4 floats and `linalg.matrix` on a real `@Matrix` of 2 to 4 by 2 to 4 take SPIR-V's own instructions by inline `asm` (design.md, "Dispatch under SPIR-V"); other sizes and complex `T` keep the inline code.

| Operation | Emitted |
| --- | --- |
| `math.dot`, `math.dot2` | OpDot |
| `math.cross` | Cross (68) |
| `length` | Length (66) |
| `distance` | Distance (67) |
| `normalize` | Length + Normalize (69), the zero test between |
| `linalg.vector.mulAdd` | Fma (50), through `@mulAdd` |
| `mix`, `mix_scalar` | FMix (46) |
| `matrix.det` | Determinant (33) |
| `matrix.inv` | Determinant + MatrixInverse (34), `error.Singular` when 0 |
| `matrix.outer` | OpOuterProduct |

- **Shaderlib's own (S2, 2026-11-12):** `reflect`, `refract`, `faceforward`, `inversesqrt`, `stepScalar`, `smoothstepScalar`, `fract` and `clamp` (on floats) are GLSL.std.450 Reflect (71), Refract (72), FaceForward (70), InverseSqrt (32), Step (48), SmoothStep (49), Fract (10) and FClamp (43) by inline `asm` under SPIR-V, the plain formulas elsewhere; `clamp` is `std.math.clamp` off SPIR-V and for integers.
- **`asm` on `@Matrix`:** the block takes the matrix type as a `"t"` input and the matrix as a `""` operand, so it never declares `OpTypeMatrix` ([#198](https://github.com/randyos-org/rstd/issues/198)).
