**Headline:** The generated reference pages (`docs/references/fitted.json`, `overview.json`, `stand-ins.json`) are rebuilt after the citation-form fix of PR 7451 and the icy-body rock fraction row of GEN.522, so the generated-files check is green again.

### Fixed

- `scripts/regenerate_generated.py --only reference-pages` was out of date on main after PRs 7450 and 7451 (a citation form in `tuning_targets.json` feeds `fitted.json`; the GEN.522 rock fraction row feeds `stand-ins.json`). The three files are the builder's own output, nothing else changes.

generation-output: unchanged
