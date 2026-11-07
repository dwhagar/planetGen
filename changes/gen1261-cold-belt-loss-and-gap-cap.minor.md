**Headline:** Belts are smaller and more honest: the solids of the cold disk lose the share that never grew to planetesimals, a gap holds no more than the planet it would have grown, and belts are made of the disk the planets grew in.

### Changed

- The solids a belt starts with (`belt_origin`, GEN.1261): past the snow line only `planetesimal_share` of the disk's solids grew to planetesimals before the pebbles drifted in (about one half for a solar disk; stand-ins BG17 to BG19), and a gap belt holds no more than the isolation mass at its centre (`isolation_cap_earth`). In a census of 240 F, G and K systems the median outer disk fell from 10 to 2 Earth masses and the largest gap belt from 3.2 to 0.6; the share of systems with a belt did not change (inner 11 percent, cold 70 percent against 14 to 20 observed).
- The debris-disk cut (`debris_dust.DETECTABLE_FRACTIONAL_LUMINOSITY`, stand-in DU1, L3) moves from 1.5e-4 to 3e-5: with the lighter cold belts only 1 to 3 percent of G, K and F dwarfs would shine above the old cut (7, 34 and 4 percent before), against the surveys' 17 percent for the three together; the new cut gives 14.7, 31 and 7 percent (300 systems each, 17.7 percent together).
- Belts are built in the system's own disk: `belt_birth.build` takes the planets' `disk_scatter` (0.7 dex) where it used the median disk, so a light disk has light belts.
- Tests for the cold share, a disk with no drift loss, the isolation cap and the disk scatter; the stand-ins BG17, BG18 and BG19 are in the GEN.725 registry.
