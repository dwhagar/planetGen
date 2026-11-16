**Headline:** A small or cold icy body's tidal heat and circularisation now use solid ice's tide instead of the ocean shell's, so frozen-through bodies are no longer heated at the shell's cap rate.

### Changed
- GEN.1759: `interior_heat.icy_regime` is now the one chooser (shell on an ocean, frozen through, or rock) that the spin (`body_spin`), the tidal heat and circularisation (`interior_heat.tidal_response`, `heat_at`, `update_heat`, `eccentricity_tide_w`) and the sea's tide (`ocean_tides`) all read. It reuses the choice 6787 added to `body_spin` (stand-ins IC6, IC7), moved so there is one copy.
- A body frozen through now gets the homogeneous ice response (`tidal_rheology.ice_response` at the mean conducting temperature) in its tidal heat and orbit circularisation. Where the shell formula gave the 1.5 cap, solid ice has 1e-4 to 1e-3: the median tidal flux of such bodies is 9.5e-10 to 4.3e-17 W/m2 on 300 seeds.
- A water-rich body whose surface is above the melting point of ice has no shell, so its sea's tide uses a rigid Love number instead of the ice shell's.
- Not changed: the tidal heat of those warm water-rich bodies (60 of 5,878 icy bodies on the census seeds) stays as before, because the orbit is circularised before the surface temperature is known (it sees a 100 K shell); that order is GEN.1760.
- Generation output changes: tidal heat, interior temperature and orbit of frozen-through icy bodies.
