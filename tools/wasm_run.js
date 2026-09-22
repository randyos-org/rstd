// Runs a WASI command module under Node, for the wasm32 test step:
//
//     node tools/wasm_run.js <module.wasm> [args...]
//
// stdio goes through, and the process exits with the module's status, so a
// failing test fails the step. Needs a Node with the `node:wasi` module.

const fs = require("node:fs");
const { WASI } = require("node:wasi");

const [file, ...args] = process.argv.slice(2);
if (!file) {
    console.error("usage: node wasm_run.js <module.wasm> [args...]");
    process.exit(2);
}

const wasi = new WASI({
    version: "preview1",
    args: [file, ...args],
    env: {},
    preopens: {},
    returnOnExit: true,
});

WebAssembly.compile(fs.readFileSync(file))
    .then((module) => WebAssembly.instantiate(module, wasi.getImportObject()))
    .then((instance) => {
        process.exitCode = wasi.start(instance);
    })
    .catch((err) => {
        console.error(err);
        process.exit(1);
    });
