**Headline:** The libraries the site ships itself are now locked: a test fails when a vendored file is edited by accident.

### Added
- **`static/vendor/VENDORED.json` (OPS.35).** For three, Shoelace, xterm and its fit addon, and the two TanStack cores: the version, the npm `dist.integrity`, the bundler, the licence and the SHA-256 of every shipped file. `scripts/vendor_lock.py` writes it (`--check` reports differences) and `scripts/vendor_shoelace.py` rewrites it after a Shoelace bump; `tests/test_vendored_lock.py` checks every file under `vendor/` against it and fails on a change or an unlisted file.

### Changed
- **The third-party notices name the bundler and the licences inside Shoelace.** esbuild 0.28.2 (it reproduces three and both TanStack bundles byte for byte), and Lit (BSD-3-Clause, Google LLC), Floating UI and `@shoelace-style/localize` (MIT), which the vendored Shoelace chunks bundle.
