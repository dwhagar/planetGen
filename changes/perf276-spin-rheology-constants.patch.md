**Headline:** A planet's spin history runs about a third faster: the rock rheology and each raiser's tide constant are worked out once for the planet, not at every one of the hundreds of thousands of torque evaluations.

### Changed
- PERF.276: `tidal_rheology.k2_over_q_function` holds the rigidity, viscosity, self-gravity pressure and Andrade constants of one body at one temperature and melt share, and gives `k2/Q` of a tide frequency by the same arithmetic in the same order as `response(...).k2_over_q` (a test holds it equal to the last bit). `spin_torques.rheology_dissipation` uses it, and `Torques` keeps `K_g` of each raiser from construction. `body_spin.evolve` over four G2V systems fell from 3.4 s to 2.3 s; the build from 43.8 s to 41.7 s. A digest of every body's mass, distance, eccentricity, day, tilt and temperature, and each moon's, is the same before and after.

generation-output: unchanged
