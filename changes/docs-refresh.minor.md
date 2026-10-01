### Changed
- **README and INSTALL split.** `README.md` now says what planetGen is,
  what it does, what it needs and how to use the website and the command
  line. The new `INSTALL.md` walks from nothing to a running site on
  Linux, Windows or macOS with the provided scripts, says where every
  example config lives, and covers updates and scheduled maintenance.
  The old `pip install .` setup step is gone: the install scripts set up
  the libraries and everything runs the checkout's code. The full
  command-line reference moved to `docs/cli.md`.
- **TODO items have permanent category IDs** (`UX.1`, `MAP.16`, ...),
  a plain running count in each category like the schema version,
  replacing the running numbers. Bugs and features are listed under the
  item they belong to, `docs/TODO.md` is the one place that links an
  item to its design document, and code tags read `TODO(MAP.16)` (a test
  checks each names an open item). `docs/design/todo-number-map.md` maps
  every old number, by date, and the short-lived dotted IDs (`MAP.2.1`)
  to the new IDs, for the changelog, commits and PRs that cite them.
- **Every reference doc checked against the code.** `database-schema.md`
  now describes schema v44 and its tables; `api.md`,
  `html-interface.md`, `config.md`, `system-file-format.md`, the
  deployment guides and the rest have their errors fixed (for example,
  `update.sh` resets to the branch tip rather than refusing to run over
  local changes).

### Added
- `docs/design/architecture.md`: how the program fits together, mapping
  every script, package and module to what it holds and tracing the main
  flows with diagrams.
- `docs/design/design-decisions.md`: the big design choices, when they
  were made, why, and what was rejected. The design documents were
  brought up to date with their reasons; superseded ones moved to
  `docs/design/archive/` and `docs/analysis/archive/`.
