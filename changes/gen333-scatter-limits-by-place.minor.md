**Headline:** The star scatter draws its line by place: the galactic core, the halo and everywhere else each have their own luminosity floor and mass limit.

### Changed
- **A band takes the core, baseline or halo tier from the term that dominates its density (GEN.333).** The nuclear cluster and disc make a band core, the halo makes it halo, everything else is baseline. The core draws at its own floor and mass limit (the baseline's times 1.1 until set: 9,900 and 15.4 at the defaults), the halo at its own floor (default 500) and the baseline's mass limit. The mass and luminosity passes, the marks that keep a sector from the luminosity pass, the single-pass redo, the staged `--bright-stars-down-to` band, a sector's backfill and mass backfill, the fill and the cell resolution all read the limits of the band they work in. A galaxy scattered before this has the wrong rows in the core and halo: `planetgen plan --redo-scatters mass luminosity` draws them again.
- **The sector page and `/sectors/<id>` name the band's tier and the limits it was scattered at.**

TODO items resolved: GEN.333. Added: none.
