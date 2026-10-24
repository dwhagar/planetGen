**Headline:** The atmosphere chemistry now has one equilibrium solver and one thermodynamic data table instead of two.

### Removed
- **The first equilibrium module (GEN.672).** `physics/equilibrium.py`, its data file `thermo_nasa7.json` and `scripts/build_thermo_table.py` held the same NASA data as `gas_equilibrium.py` with `species_thermo.py`. The three hot-Jupiter regime edges (methane to carbon monoxide, water survival, hydrogen dissociation) now run on the one solver in `physics/equilibrium_edges.py` and give the same temperatures; no generated value changes.
