**Headline:** Every galaxy's central black hole and nuclear star cluster come from its own mass, bulge and seed, not from fixed Milky Way figures.

### Changed
- **The central black hole follows the galaxy** (`galaxy/structure.py`, `galaxy/density.py`). A galaxy with a classical bulge (the early stage) takes the bulge relation, `M_BH / 1e9 = 0.49 (M_bulge / 1e11)^1.16` with 0.29 dex of scatter (Kormendy and Ho 2013); any other galaxy takes `log M_BH = 7.45 + 1.05 log(M* / 1e11)` with 0.24 dex (Reines and Volonteri 2015). Each galaxy's draw about the relation comes from its seed. The bulge's velocity dispersion follows from the mass by the M-sigma relation (Kormendy and Ho 2013, Eq. 3). A black hole or a quasar made with no galaxy takes the hole of a randomly drawn host galaxy.
- **The nuclear star cluster follows the galaxy** (`galaxy/structure_draw.py`, GX3). Its mass is `log M = 6.51 + s log(M* / 1e9)` with 0.6 dex of scatter (Neumayer, Seth and Boeker 2020), and the slope `s` is drawn per galaxy uniform on 0.48 to 1.36, the Neumayer and Georgiev values both kept (Boss). The share of galaxies with a cluster is a labelled placeholder (0.85 up to 1e9 solar masses, falling to 0.2 at 1e11) until the occupation fractions are read; a galaxy without one has no nuclear disc either. The nuclear disc and size ratios stay Milky Way stand-ins.
- Removed from `tuning.py`: `NUCLEAR_CLUSTER_BULGE_RATIO`, the supermassive black hole mass range and dispersion, and the quasar mass range. The structure gains `stellar_mass_solar`, `classical_bulge`, `black_hole_deviation` and `nuclear_cluster_mass_solar`.
- The generator epoch is 25.
generation-output: changed

### Added
- `tests/test_central_black_hole_and_cluster.py` (10 tests): determinism, both black hole relations, the mass trend over 1,200 seeds, M-sigma, the cluster's power law and slope range, the nucleation share, the density terms without a cluster, a hole made with no galaxy, and that the old constants are gone.

TODO items resolved: GEN.402. Added: none.
