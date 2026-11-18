**Headline:** One routine builds the planets of every host star, instead of four copies of the same seven steps.

### Changed
- `StarSystem._build_host` draws, spaces, trims, adds belts and star history to, and widens for mass loss the planets of one host (GEN.1656). The primary, a wide pair's secondary, a tree member and a node's proxy each had their own copy of the sequence; each fix (the ceiling trim, the white dwarf widening, the belts) had to be made four times. The output is unchanged: 42 systems (single, wide, member and node hosts) hash the same before and after.

generation-output: unchanged
