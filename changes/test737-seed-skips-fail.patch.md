**Headline:** Tests that skipped when their fixed seeds drew nothing now fail, so a draw that moves cannot hide a test; the two that were skipping today (white dwarf sky text, moon ocean search) run.

### Changed
- **TEST.737: 34 "no X in these seeds, draws or sample" skips become failures, and the two that were skipping on the current draws now run.** A survey of the 34 files that hold these skips found only two firing: the white-dwarf sky-text test (it now asks for planets, so the draw always has some) and the moon ocean-depth range search (it needs three values, which the sample galaxy holds). The conditional skips for a missing browser feature (WebGL), a missing data set, and an out-of-range shell and band stay. Tests only; no generated output moves.

generation-affects: none
