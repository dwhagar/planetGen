**Headline:** The wiring register now says why the Hycean box and the shoreline module are not read by the generator, and which item reads them.

### Changed

- GEN.1523: `hycean_candidate` and its row HC1 are marked by design (W26 has no occurrence rate; the `hycean` state follows from the generated planet, WR17; the box is for the census and tests). `atmosphere_retention` and rows SH1, SH2 and SH4 now name GEN.1548 (the air stage rebuilt once on the final orbit reads the shoreline there); SH3 names GEN.1558 (the caps are retired there); the energy-limited term is not wired from the module (GEN.1558 reuses the thermosphere code). `small_body_composition.draw_ring` was wired by `rings.build_ring` (part 1, PR 5990); `helium_fractionation` is GEN.1527. No code changes.

generation-output: unchanged
