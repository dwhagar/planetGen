**Headline:** A star near the hydrogen-burning limit no longer draws a companion lighter than any brown dwarf can be.

### Fixed
- **A multi-star tree failed for a primary near 0.08 solar masses.** The other members took `m1 * U(0.1, 1)`, which reaches
  0.0085 solar masses for a primary of 0.085, below the lowest cooling track (0.01) and the deuterium limit of 13 Jupiter
  masses (`brown_dwarf.mass_range_solar`), so `brown_dwarf.state` raised "no cooling track". The mass ratio is now no lower
  than puts the member at the deuterium limit; lighter objects are planets, filled in as planets are (GEN.596). A primary
  of 0.124 solar masses or more draws as it did.

generation-output: changed
