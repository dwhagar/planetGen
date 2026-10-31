**Headline:** The mantle viscosity draw uses the laboratory's activation energies (260 to 390 kJ/mol), so the temperature dependence of a planet's tidal heating is steeper than before.

### Changed
- **GEN.1185.** `INTERIOR_ACTIVATION_RANGE_KJ_MOL` goes from 200 to 300 to 260 to 390 kJ/mol (diffusion creep, Korenaga and Karato 2008 Table 3) and `INTERIOR_VISCOSITY_PREFACTOR_RANGE_PA_S` from 2.6e10 to 5.3e13 to 9.3e7 to 1.9e11 Pa s, the old bounds times 3.56e-3, so the viscosity of a mantle at 1,600 K in the middle of the ranges holds. The same draws give each body a different `mantle_eta0_pa_s` and `mantle_activation_kj_mol`, so interior heat, tidal heating, spin-down and anything that reads the mantle's viscosity change a little: bodies far from 1,600 K change most. Stand-in rows IHA1 (prefactor, tier T4, tweak L3) and IHA2 (activation, T1, L1) added in section 6.163. Io's tidal power over the corners stays 9e12 to 5e13 W.

generation-output: changed

TODO items advanced: GEN.1185.
