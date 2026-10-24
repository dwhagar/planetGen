**Headline:** Every stand-in value now has a source tier and a tweak level, and one generated page lists them all.

### Added
- **Stand-in tiers and tweak levels (GEN.725).** The stand-in registry of the gap-filling plan carries two new columns on every row: the source tier (T1 verified down to T9 educated guess based on needed output) and the tweak level (L0 as envisioned up to L4 fitted to outcome). `docs/design/research/stand-in-values.md` lists every row, built from the registry by `scripts/build_stand_in_values.py`, and a test fails when a row lacks either column or the page differs from the registry. No generated value changes.
