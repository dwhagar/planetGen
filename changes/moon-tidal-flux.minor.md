**Headline:** The tide that heats an icy moon now feeds the ocean's life under its ice, so Europa- and Enceladus-like moons score on their real heat.

### Changed
- **Tides in the chemical energy of an ice-covered ocean.** `hydrosphere.generate_hydrosphere` kept the tide of an icy moon only long enough to thicken or thin its ice lid and then threw it away; the chemical energy term of the habitability score read the radiogenic heat alone, one to three orders of magnitude too low for a Europa or an Enceladus. The tide (`tidal_flux_w_m2`) is now kept on the moon, and `habitability_world.generate` adds it to the moon's own heat flow before the energy term (1 percent of the heat flow, as for any body).
- **Migration 0152, schema v152.** `planets` and `moons` gain `tidal_flux_w_m2` (DOUBLE, nullable; a planet keeps NULL). The forced eccentricity behind a tide is not stored, so a moon saved earlier keeps NULL and its energy flux; only a moon generated now has it. Generator epoch bumped.
- **A moved planet's moons are scored again (GEN.704).** `validation.reconcile_moved_planet` refreshed a moon's air for its parent's new distance from the star but kept the dose, UV and habitability it had drawn at the old one; each moon that keeps its class now gets `generate_spin`, which ends with them, as a moved planet does.
- **The property sheet shows the tide.** `describe/sheet.py` lists `tidal_flux_w_m2` in the Water group of a planet and a moon.

generation-output: changed
