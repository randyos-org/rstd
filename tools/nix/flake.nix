{
  # QEMU for the qemu rigs on another machine (`-Dtest-host`), pinned by
  # flake.lock; `nix flake update` moves it on. tools/test_host.zig copies
  # this directory to the host and runs run.sh in its shell.
  description = "rstd's qemu rigs' QEMU";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";

  outputs =
    { nixpkgs, ... }:
    {
      packages = nixpkgs.lib.genAttrs [ "aarch64-darwin" "x86_64-darwin" "aarch64-linux" "x86_64-linux" ] (
        system:
        let
          pkgs = nixpkgs.legacyPackages.${system};
        in
        {
          # every system target, with its firmware
          inherit (pkgs) qemu;
          # run.sh's slots
          inherit (pkgs) flock;
          default = pkgs.qemu;
        }
      );
    };
}
