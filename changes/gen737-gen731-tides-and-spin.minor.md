**Headline:** The tide of a body's real orbit heats it and wears down its spin: Io-like moons melt and radiate, Mercury-like worlds are caught in a resonance, hot Jupiters lock, and the Moon-like ones stay cold.

### Changed
- **Tidal heating comes from the orbit's real eccentricity (GEN.737).** A rocky or icy planet or moon is heated by the eccentricity tide
  of the orbit stored with it (Peale, Cassen and Reynolds 1979), from the eccentricity of GEN.727. The drawn forced eccentricity and the
  random 1 to 50 times resonance boost of small icy moons are gone, and so is the tide of an icy moon only: planets and rocky bodies have
  one too.
- **A rocky body's rheology has two states.** Solid, its Love number comes from its rigidity (the Moon's 0.021) and Q is 100, which gives
  the Moon 2e8 W. If that heat is at least the body's own and would melt its mantle, the melt erupts and carries the heat out, the mantle
  holds the melt temperature and the partly molten body is flexed far more: Io comes out at 9e13 W, 2.2 W/m2.
- **The heat is stored and counted.** The tide goes into the interior heat model (a heat source of the mantle, and the melt state), is
  stored in `tidal_flux_w_m2` for every rocky and icy planet and moon, and is added to the body's own heat flow in the habitability
  energy term and the ice-shell thickness. The own heat flow (`internal_heat_flux_w_m2`) and the Urey ratio leave the tide out.
- **Enceladus makes the tide alone.** Its equilibrium tide is about 0.4 GW; the 16 GW it radiates is an open problem in the literature
  that the old boost hid.
- **A planet's day is its birth spin worn down by tides (GEN.731).** A rocky planet is born spinning at the rate and tilt the giant
  impacts that built it leave (Kokubo and Genda 2010: Gaussian components, a median day of 2 to 3 hours, an isotropic axis), a moon or a
  giant from its own range, and the tide of what it orbits then slows the spin for the age of the star, with the body's own Love number and
  quality factor (the same that circularise its orbit). A close body ends locked, one turn per orbit with its axis upright; an eccentric
  one is caught in a spin-orbit resonance (Mercury's 3:2 at its eccentricity of 0.2; windows from Makarov 2012); a far one keeps its spin.
  A hot Jupiter locks in about a hundred thousand years; the Moon, in tens of millions. This replaces the flat 8 to 48 hour day of rocky
  planets, so typical planetary days are now much shorter than the Solar System's (24 hours): the source simulates the last giant
  impacts only.
- **Stand-ins recorded.** Rows TH1 to TH9 and SP1 to SP10 in the stand-in registry (GEN.725).
generation-output: changed
