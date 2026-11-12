**Headline:** The unused dwarf-population module is retired: what it did, the belt member status words and the size-law order statistics, is already done by the running code.

### Changed
- GEN.1092 (wiring triage, GEN.1559): `physics/dwarf_population.py` (GEN.1090) was built and read by nothing. Its status word, Soter's mu, Margot's Pi and round mass are `physics/body_status.py`, which the belt page already uses for every stored member (the registry rows BR6 and BR7 now name `body_status`), and its counts and order statistic are the belt ledger's own draw (`physics/belt_ledger.py`). The module and its test are removed, with its line in the wiring register. Class D for a belt's heaviest round bodies comes from the planet generator when GEN.1510 makes them planets, not from a second rule.
- No generated value changes.

generation-output: unchanged
generation-affects: none
TODO items resolved: none
