**Headline:** The golden galaxy now has a digest for each generation stage, and a train can check which stages its output moved against the stages its notes declared.

### Added
- **`tests/golden_stages.py` and `scripts/train.py stages`** (OPS.147): the golden test writes one digest per stage (plan's own rows, scatter tables, fill sections, facilities) to the file named by `PLANETGEN_STAGE_DIGESTS_OUT` when it builds the pinned galaxy, with no second build. `train.py stages --before A --after B` fails a stage that moved when no declared stage covers it (a declared stage covers itself and the stages after it) and warns about a declared stage that did not move. The golden files and the epoch are untouched. Documented in `docs/reference/testing.md`.

generation-output: unchanged
generation-affects: none
