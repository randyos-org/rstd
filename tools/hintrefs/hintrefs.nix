# hintrefs.c against nixpkgs' FreeType (tools/nix/hintrefs.sh)
{
  runCommandCC,
  freetype,
  pkg-config,
}:
runCommandCC "rstd-hintrefs"
  {
    nativeBuildInputs = [ pkg-config ];
    buildInputs = [ freetype ];
  }
  ''
    $CC -O2 -std=c11 -Wall -Wextra -Werror ${./hintrefs.c} -o $out $(pkg-config --cflags --libs freetype2)
  ''
