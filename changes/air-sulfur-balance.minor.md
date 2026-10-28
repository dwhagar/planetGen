**Headline:** The sulfur a world keeps in its air is what the sulfate and sulfide sinks of its crust leave, instead of a drawn share.

### Changed
- **The sulfur balance replaces the stand-in (GEN.815).** `physics/surface_air.py` gives the mixture draw the sulfur that `physics/sulfur_balance.py` leaves in the air: the sulfur the interior degassed, less the sulfate (calcium) and sulfide (iron) the crust it built can bind, at the removal lifetime of rain-out over an ocean or of sulfation of the rock on a dry surface, split by the redox of the gas. Sulfur beyond what the crust holds stays in the air whole. A wet world no longer holds hundreds of bar of sulfur gas, and a dry one keeps almost none unless its inventory outgrows its crust.

### Removed
- **`SULFUR_SURFACE_KEPT` (1e-5 to 1e-2 of the degassed sulfur kept with an ocean)**, its drawn constant `sulfur_surface_kept` and its registry row (AM80).

generation-output: changed
