**Headline:** An ocean now holds some of its world's nitrogen dissolved: the air solve takes the nitrogen out by Henry's law, as it already did the sulfur gases.

### Added
- GEN.1834: `mixture_sinks.remaining_fraction` (renamed from `sulfur_remaining_fraction`, which had no other caller) takes `N2` beside `SO2` and `H2S`, and `surface_air._dissolve` (renamed from `_dissolve_sulfur`) takes it out of the column over a surface ocean. The law is the verified N2 entry of the gas table (kH 0.000625 mol/(kg bar) at 298 K, 273.15 to 323.15 K, the nearest end outside), over the same ocean column the sulfur uses (all the surface water); no new constant. The share does not depend on the pressure: an Earth ocean holds 0.6 percent of the air's nitrogen, one ten times as deep 6 percent. The nitrogen that crust, sediment, life and subduction keep is not added: it is inside the degassing efficiency NB3 (GEN.1833, `reports/gen1833-nitrogen-sinks-2026-11-22.md`). Nitrogen takes nothing outside the law's 273.15 to 323.15 K (a frozen sea is cut off from the air, and the law carried to a hot or frozen water world would take its whole air); the sulfur gases keep the nearest-end clamp. The salt of the sea (about 7 percent less solubility) is left out.

generation-output: changed
generation-affects: fill
