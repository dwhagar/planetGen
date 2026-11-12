**Headline:** The one-off system page and `planetgen system` can set the weirdness dial and force particular oddities.

### Changed
- **Weirdness dial and Force oddities (UX.299).** `planetgen system` takes `--weirdness-dial X` (0 turns every oddity off, 1 is the table's own rates) and `--force-oddities W10,G1` (the ids of `physics/oddities.py`; a body whose conditions allow them gets them without the roll, and the gate still applies). Both set the existing `SystemConfig.WEIRDNESS_DIAL` and `FORCE_ODDITIES` (GEN.1047); an unknown id or a negative dial is refused. The one-off system page has the same two boxes under Options, with the ids listed. Nothing stored changes.

generation-output: unchanged
