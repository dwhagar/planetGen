**Headline:** Every star now draws its birth [Fe/H] first, and its lifetime, subgiant and giant branch and horizontal-branch colour follow that composition.

### Changed

- **A star's composition comes before the star (GEN.651).** `birth_fe_h` draws a star's [Fe/H] from its population at its place (or its cluster's or satellite's own distribution) ahead of every other draw. The main-sequence lifetime (Hurley, Pols and Tout 2000), the end of the subgiant phase, the red giant branch and its tip, the Wolf-Rayet threshold and the horizontal-branch colour all read it, so a metal-poor star lives and brightens differently from a metal-rich one of the same mass. `star_params` carries it as `metallicity_fe_h`.
- **The bright-star sampler keeps one table per composition.** The cell bounds of a sampling table hold only for the composition they were built at, so `star_population` builds a table for each 0.25 dex bin of [Fe/H] (`tuning.BRIGHT_TABLE_FE_H_STEP_DEX`), takes each cell's bound over the bin's two edges and its centre with the stretch of every phase window (`tuning.BRIGHT_TABLE_BOUND_MARGIN`), and draws a star's composition in proportion to the stars each composition can make (`_draw_fe_h`), so a metal-poor place no longer fills its bright giants with stars its composition cannot make. The workers warm every bin.
- **Stand-ins until GEN.661 and GEN.662.** The expected count of bright stars still uses the share at the reference composition, and a cluster or satellite star is drawn at the reference composition rather than its component's value. The lifetime fit's metallicity is held to the range the Hurley, Pols and Tout formulae were fitted over (Z 0.0001 to 0.03), so a very metal-poor star no longer extrapolates.

generation-output: changed
