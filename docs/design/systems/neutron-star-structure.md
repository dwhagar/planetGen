# Neutron star structure: a TOV integration over the 21 equation-of-state tables (GEN.968)

Status: lock-free step built (GEN.968, parent GEN.743). Nothing in generation calls it yet; the wiring is GEN.969 and rides the star-epoch merge train. Code: `physics/neutron_star_structure.py`, `scripts/build_neutron_star_table.py`. Data: `datasets/neutron_star_matter.py` (Malik, Pais and Providencia 2024, 21 unified equations of state, GEN.548) and the reduced `physics/data/neutron_star_curves.json`. Tests: `tests/test_neutron_star_structure.py` (synthetic polytropes and the committed table against NICER and the pulsar masses; no download) and `tests/test_neutron_star_structure_real.py` (the table against a fresh integration of the real tables).

## What it does

A neutron star is fixed by its equation of state and one number, the pressure at its centre. The module integrates the Tolman-Oppenheimer-Volkoff equations in the logarithm of the pressure, from the centre to the lowest pressure of the table (about 1e4 g/cm^3, the surface), and answers:

- `structure(eos, central_pressure)`: mass, radius, compactness, the radius where the crust begins and where neutron drip is, and the mass of the crust.
- `mass_radius_curve(eos)`, `maximum_mass_solar(eos)`: the stable branch of one equation of state and the heaviest star it holds.
- `radius_km(eos, mass)`, `supporting(mass)`, `radius_band_km(mass)`: the radius of a star of a given mass, the equations of state that hold such a star, and the lowest, median and highest radius across them.

## What the 21 tables give (as built, 10-30)

| Quantity | The 21 tables | Measured |
|---|---|---|
| Maximum mass | 2.01 to 2.75 solar masses (19 of 21 above 2.08) | PSR J0740+6620, 2.08 +- 0.07 (Fonseca et al. 2021) |
| R(1.4) | 11.8 to 13.8 km, median 12.9 | 12.45 +- 0.65 (Miller et al. 2021) |
| R(1.34) | median 12.9 | 12.71 +1.14 -1.19 (J0030+0451, Riley et al. 2019) |
| R(2.08) | 11.5 to 14.1 km, median 12.4 | 12.35 +- 0.75 (Miller et al. 2021); 13.7 +2.6 -1.5 (J0740) |
| Crust of a 1.4 solar mass star | 0.4 to 2 km thick | about 1 km |

The stiffest tables (EOS 1 to 5) put R(1.4) above the Miller band by up to two of its widths. GEN.969 should draw a star's equation of state with a weight from these constraints (NICER radii, the 2.08 solar mass floor) and not uniformly; the module gives the radius of each table and leaves the weight to the wiring.

## Numerical choices

- Units are G = c = 1 in kilometres; 1 MeV fm^-3 is 1.3238e-6 km^-2 and a solar mass is 1.4766 km.
- The equation of state between rows is PCHIP in log e against log P (no overshoot), resampled on 6,000 points for speed.
- DOP853 with a relative tolerance of 1e-6; against 1e-10 the radius and mass of a star change by about 1e-5.
- Checked against the Newtonian n = 1 polytrope (R = sqrt(pi K / 2G) within 0.5%) and the Buchdahl bound (compactness under 4/9).
- Cost: about 1 s a table for its maximum mass and 50-point curve, about 22 s for all 21. Generation does not integrate: `scripts/build_neutron_star_table.py` writes `physics/data/neutron_star_curves.json` (25 stars of each table, 28 kB: mass, radius, crust thickness and crust mass), `radius_km`, `crust_thickness_km`, `maximum_mass_solar`, `supporting` and `radius_band_km` read that, and the real-data tests hold it to the integration. `--check` fails if the file is out of date; rebuild it when the integrator or the crust density changes.

## Stand-ins (registry 6.32, GEN.725)

The crust-core density 0.08 fm^-3 and the neutron-drip density 2.6e-4 fm^-3 mark the layers; the tables do not. NS1 to NS3.

## Open for GEN.969 (the wiring)

- `physics.neutron_star.radius_km(mass, deviation)` keeps its signature: the standard-normal `deviation` the caller draws (one `draw.gauss()` per star, `generation/phenomena/compact_remnant.py`) becomes a quantile that picks the equation of state, so no draw is added or removed. The 21 are weighted by how well each fits NICER (a normal weight on R(1.4) about 12.45 +- 0.65 km, Miller et al. 2021, and the 2.08 solar mass floor of J0740) and by whether it supports the star's mass (`supporting`). That is a stand-in to register (tier T2, tweak L1) when it is wired.
- A star's radius comes from `radius_km(eos, mass)`; its crust (for the surface and the thermal model) from `crust_thickness_km`. `tuning.NEUTRON_STAR_RADIUS_ANCHORS` and `NEUTRON_STAR_RADIUS_CLIP_SIGMA` go (no backward-compatibility shim).
- The mass range of a neutron star (`tuning.NEUTRON_STAR_MASS_RANGE_SOLAR`, 1.1 to 2.2) stays; the heaviest equations of state allow up to 2.75.
- This moves every neutron star's radius, so it is `--affects fill` and rides train 3.

The helium-core and oxygen-neon white dwarf layers (the other half of GEN.743) are not in this step; they need a white dwarf structure source the Bedard et al. sequences do not hold.
