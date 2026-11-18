**Headline:** the greenhouse table's fit script reads the table the way the generator does, and gains the cold dense families that show where the table is weak; the table file is unchanged.

### Changed
- GEN.1796 (tools only): `scripts/fit_greenhouse_table.py` `validate` now reads each air through `greenhouse_table.warming` (minor gases counted by share, the haze guard, the thin-air and deep-air rules) instead of the bare network, so the stated accuracy is what the product has (the first fit stated it so; round 3 of GEN.1780 had read the bare network, 79 K against 55 K for the same file). New: `cold --cia` and `cold --gen` (cold dense nitrogen airs of light worlds and with the generator's own minor gases), `nullgas` (census airs holding a gas with no ir flag, held at the starting file's answer), `solve --cold-dense` (the held airs the starting file had wrong), `fit --boost`. Findings and the round files are in `/mnt/project-files/reports/gen1796-refit-2026-11-18/`; the table still misses the 50 K clear-air bound on its stated set and the refit is open under GEN.1796 and GEN.1781.

generation-output: unchanged
generation-affects: none
