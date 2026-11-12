**Headline:** A version can no longer be stamped while built physics sits unread: the wiring guard has a release mode that fails on any unread item the register does not mark by design.

### Changed

- TEST.718: `python scripts/wiring_guard.py --release` prints every module, function, registry row, data file, dataset and column that is built and read by nothing and that `docs/design/research/wiring-register.md` does not mark `by design`, and exits 1 while there is one (124 today; the burn-down is GEN.1559). It is the third release gate beside the SEC.33 and TEST.261 audits and is documented in `docs/plan/execution-notes.md`. Two tests pin it: a seeded unread item fails the gate and an empty list passes.

generation-output: unchanged
generation-affects: none
