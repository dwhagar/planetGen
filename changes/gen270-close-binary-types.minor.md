**Headline:** Close binaries now come in their observed kinds: contact binaries, Algols, blue stragglers, merger products, Be stars with a stripped companion, and the rare R Coronae Borealis, FK Comae and Thorne-Zytkow stars.

### Added
- **Contact binaries and Algols follow their gates.** A contact binary of a solar-type primary needs an age of at least 0.1 Gy, an Algol an age under 3 Gy, and neither holds a contracting or substellar object. The donor of an Algol that was the heavier at birth carries `mass_reversed`, and the lighter star is a `mass_gainer`. The tight tail of the period law is calibrated (`tight_tail_weight` 6.0): about 0.5% of FGK systems hold a contact binary, and 60% of contact binaries have a tertiary.
- **Merger products and mass gainers.** 8% of young massive stars are merger products with one star fewer and a rejuvenated age, 0.5% of older stars are blue stragglers, and 20% of young massive stars are mass gainers.
- **Be stars.** 20% of young B stars are Be stars, and 30% of those have a stripped helium companion that spun them up.
- **R Coronae Borealis, FK Comae and Thorne-Zytkow stars.** Old giants take these at the catalogue's Galaxy counts; a Thorne-Zytkow object is a single member flagged theoretical.
- **Eclipsing and spectroscopic pairs** are derived from a tree (`MemberTree.eclipsing_nodes`, `spectroscopic_nodes`), not stored.

### Changed
- **Generated systems change.** The period law, the gates and the new draws move the draws of star systems, so every galaxy regenerates differently (generator epoch recorded).

generation-output: changed
