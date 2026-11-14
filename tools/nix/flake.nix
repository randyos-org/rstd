{
  # QEMU for the qemu rigs on another machine (`-Dtest-host`), pinned by
  # flake.lock; `nix flake update` moves it on. tools/rig/test_host.zig copies
  # this directory to the host and runs run.sh in its shell.
  description = "rstd's qemu rigs' QEMU";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  # unstable dropped x86_64-darwin (26.11); 26.05 is its last release,
  # supported to the end of 2026
  inputs.nixpkgs-x86_64-darwin.url = "github:NixOS/nixpkgs/nixpkgs-26.05-darwin";

  outputs =
    { nixpkgs, nixpkgs-x86_64-darwin, ... }:
    {
      packages = nixpkgs.lib.genAttrs [ "aarch64-darwin" "x86_64-darwin" "aarch64-linux" "x86_64-linux" ] (
        system:
        let
          pkgs = (if system == "x86_64-darwin" then nixpkgs-x86_64-darwin else nixpkgs).legacyPackages.${system};
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
