# Eigen behind blas_bench's Fortran names (eigen_shim.cpp), a benchmark
# baseline only: compiled against nixpkgs' upstream Eigen by nixpkgs' Zig, for
# Zig target `target` and CPU `mcpu`, on whatever system builds it, so the
# Mac builds the Linux ones too. tools/nix/blas.sh builds it.
{
  stdenvNoCC,
  zig,
  eigen,
  target,
  mcpu,
}:
let
  ext = if builtins.match ".*-macos.*" target != null then "dylib" else "so";
in
stdenvNoCC.mkDerivation {
  pname = "rstd-bench-eigen-${target}";
  inherit (eigen) version;
  dontUnpack = true;
  nativeBuildInputs = [ zig ];
  buildPhase = ''
    runHook preBuild
    export ZIG_GLOBAL_CACHE_DIR=$TMPDIR/zig-global ZIG_LOCAL_CACHE_DIR=$TMPDIR/zig-local
    zig c++ -target ${target} -mcpu=${mcpu} -O3 -DNDEBUG -std=c++17 -shared -fPIC -fvisibility=hidden -s \
      -I${eigen}/include/eigen3 ${./eigen_shim.cpp} -o libeigen_shim.${ext}
    runHook postBuild
  '';
  installPhase = ''
    runHook preInstall
    install -Dm755 libeigen_shim.${ext} $out/lib/libeigen_shim.${ext}
    runHook postInstall
  '';
}
