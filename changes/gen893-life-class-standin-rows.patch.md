**Headline:** Every number in the life class rules now has a row in the list of stand-in values, with its source tier and tweak level.

### Changed
- **Stand-in registry rows for the life class rules (GEN.893).** Section 6.19 of the gap-filling plan has one row for each of the fifteen life classes (solvent factor, energy band, gates, windows, supply placeholder, with the tier and tweak of the least certain number) and rows for the complexity thresholds, the hydrogen floor of non-water classes and the reference star. The generated list of stand-in values (`stand-in-values.md`) picks them up, and a new test fails when a class has no row or a row claims a better tier than its numbers.
