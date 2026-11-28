**Headline:** The reduced sulfur species now need the star's real ultraviolet-C flux, so quiet stars no longer pass the gate for free.

### Fixed

- GEN.1872: `sulfur_cloud_set` read the star's ultraviolet-C flux at the body as the flux it receives times `uv_c_share`. Since GEN.971 `uv_c_share` is the star's UV-C power as a multiple of its blackbody's (values up to about 9,800 for a quiet M dwarf's chromosphere), not a share of the light, so gate 3 (`uv_c_min_w_m2`, 0.5 W/m2) passed every star and V4 and V5 landed on bodies their own basis excludes. The flux is now what `body_radiation` already gives the surface spectrum: the bolometric flux times the blackbody fraction of 100 to 280 nm at the star's temperature times the star's UV-C multiple, through `stellar_spectrum.irradiance_w_m2` with the same library shares (`sulfur_cloud_set.uv_c_flux_w_m2(star, flux)`; reused, not rewritten). `chemistry` no longer reads a `uv_c_share` out of the values; it takes the flux from its caller or leaves gate 3 shut.
- Census of the gate on 30 stars of each type at the warm edge of their habitable zone (old reading, new reading): M5V 30/30 and 25/30 pass; M2V 30/30 and 14/30; K5V 30/30 and 5/30; K2V, G2V, F7V and A5V pass 30/30 under both. At a flux of 2,600 W/m2 (a Venus-like body round any star): M5V 30 and 29, M2V 30 and 21, K5V 30 and 8, the rest 30 and 30. Median UV-C flux at the warm edge: M5V 3.2, M2V 0.47, K5V 0.10, K2V 1.3, G2V 9.6, F7V 24, A5V 169 W/m2.
- `uv_c_min_w_m2` (0.5 W/m2) is kept at tier T9, tweak L3; the number was set on the old reading, and the new reading puts the cases its basis names where it says (a young active star at Venus's distance passes with some watts; a quiet star at its warm edge does not, with a tenth of a watt or less). Its registry basis and the CG2 row now state the new definition.

generation-output: changed
generation-affects: fill
