# The atmosphere data model

What a planet or moon stores about its air, which code makes each value, and where the numbers behind it come from. The
code is the authority; every module named here opens with a docstring that gives its sources and limits, and the designs
are in [docs/design/composition/](../design/composition/) (the colour chain's plan is
[atmosphere-colour.md](../design/composition/atmosphere-colour.md)). The colour itself is described in
[atmosphere-colour.md](atmosphere-colour.md); the columns' types and versions are in
[database-schema.md](database-schema.md). Planets and moons store the same values in the same columns; a moon's starlight is
its planet's.

## 1. The order of the work

`physics/planets.calculate_atmospheric_conditions` makes the air of a body in this order:

1. **Albedo.** The Bond albedo is drawn from the range of the body's regime (`climate_draws.albedo_range`), and the temperature
   with no air is `T_e = ((1 - A) S / (4 sigma))^(1/4)`, never below the 2.725 K background.
2. **The air and the temperature together** (`physics/surface_air.solve`, designed in
   [surface-pressure-supply-law.md](../design/composition/surface-pressure-supply-law.md)): the mixture draw, the pressure the
   interior has degassed less the sinks, the freeze-out of what the temperature condenses, the density and scale height of what is
   left, and the greenhouse that gives the temperature, repeated until the pressure and the temperature stop moving.
3. **Storing the air** (`physics/atmosphere.store_air`): the partial pressures, the mantle redox and the words of `atmosphere`.
4. **Water and ground** (`hydrosphere`, `surface_composition`), then **cloud, haze and dust** (`physics/aerosols.generate_aerosols`).
5. Later in the body's generation: **the ozone column and the doses** (`physics/radiation.generate_dose`) and **the colour**
   (`physics/body_colour`).

A body with no air, or whose air is lost, has the temperature the starlight alone gives, no gas rows, and `NULL` density,
molar mass and scale height.

## 2. What is stored

| Column | Unit | What it is | Made by |
|---|---|---|---|
| `planet_gases` / `moon_gases` rows (`gas`, `partial_pa`) | Pa | One row for each gas above zero, its formula in the gas table; the balance gas of a class (helium, sodium, ammonia) is a row like the others. No rows: no air, or a row from before v155 | `surface_air.solve`, kept by `atmosphere.store_air`; read through `mixture.partial_pressures_pa` |
| `atmospheric_pressure_pa` | Pa | The surface pressure: the sum of the partial pressures, after the freeze-out | `surface_air` |
| `atm_molar_density` | kg/mol | Despite the name, the **mean molar mass** of the mixture, `M = sum x_i M_i` | `atmosphere_state.air_state` |
| `atm_density` | kg/m3 | Surface density `rho = P M / (R T)` (ideal gas) | `atmosphere_state.air_state` |
| `scale_height_km` | km | `H = R T / (M g)`, so that `P = rho g H` | `atmosphere_state.air_state` |
| `surface_temperature_k` | K | `T_s = T_e (1 + G)^(1/4)` for the stored albedo and greenhouse quantity | `surface_air` (the greenhouse is `atmosphere_state.greenhouse_warming`, Robinson and Catling, on the infrared opacity of `physics/infrared_opacity.py`) |
| `bond_albedo` | fraction | The Bond albedo the temperature is made of | `climate_draws` (GEN.412) |
| `greenhouse_factor` | none | The quantity `G = (T_s / T_e)^4 - 1`, 0 for no air | `climate_draws.record` |
| `albedo_basis` | index | Where the two came from: 0 recorded at generation, 1 the exact inverse of a stored row, 2 a seeded draw from the class range | `climate_draws` (the migration of v142 filled old rows) |
| `mantle_redox`, `mantle_delta_iw` | class, dex | The oxygen fugacity of the mantle that outgassed the air, log10 units above the iron-wuestite buffer ([R1951] scale); `NULL` for a giant | the mixture draw (`MixtureDraw.redox_diw`, recipe rows a1 to a4, mass- and age-dependent, [R4154]) |
| `atmosphere` | text | The air in words, from the partial pressures | `atmosphere.describe` |
| `ozone_column_du` | Dobson units | The ozone column of a body with oxygen | `ozone.draw_column_du`, called by `radiation.generate_dose` |
| `ozone_loss_flag` | flag | True when there is ozone and the star's lethal event rate or the solar particle dose through the ozone layer is above its limit | `radiation.generate_dose` |
| `tau_haze`, `tau_dust`, `tau_cloud` | optical depth at 550 nm | Haze, dust and in-cloud depth (the cloud is grey) | `aerosols.generate_aerosols` (DB.68) |
| `cloud_cover` | fraction | The share of the sky under cloud, 0 to 1 | `aerosols` |
| `aerosol_draw` | 0 to 1 | One number that picks every within-bracket optical property (dust radius and asymmetry, haze exponent), so the five values repeat the optics with no other column | `aerosols` |
| `atmosphere_colour` | rgb24 | The colour of the air as `red * 65536 + green * 256 + blue`, `NULL` with no air or no star temperature | `body_colour` (DB.96) |

A column that a row made before its version does not have is `NULL`, and the page shows nothing for it.

## 3. The gas table

`physics/gases.py` is the one constants table of the gases and volatile compounds the chain reads: 112 rows, loaded from
`physics/data/gases.json`, which `scripts/build_gas_table.py` builds from the tables of the research parts
([atmosphere-composition-parts](../design/composition/atmosphere-composition-parts/), parts 01 to 03). A row holds the
formula (the key everything cites), name, CAS number, molar mass (kg/mol, summed from the element table; the research's own figure
is kept for the audit), the infrared class (`ir`, `cia`, `weak`, `none`), the volatility class, the lowest photolysis threshold and
its channel, whether it absorbs near ultraviolet, triple, boiling and critical points, ideal-gas heat capacity at 100, 300 and
1000 K, the latent heats, and where it has them the Rayleigh cross sections at 450, 550 and 650 nm, the Henry constant, the
Lennard-Jones constants that give the viscosity, the infrared bands and the ultraviolet bins.

Every row says where it comes from (`source`) and how far it is checked (`state`: `verified`, `partly verified` with its `limit`,
or `open`). The loader refuses a file or row without those fields, a duplicate formula, a formula it cannot read, and a molar mass
of the research part that is more than 0.01 g/mol from the element sum. It has no weights and no class gate: the recipes that
weight the gases (`physics/atmosphere_recipes.py`) cite its rows by formula. Units: kg/mol, K, Pa, J/(mol K), kJ/mol, cm2 per
molecule, mol/(kg bar), uPa s. A number the research does not give is `None`.

## 4. The mixture and how the values agree

The mixture is a plain dict of gas formula to partial pressure in Pa (`physics/mixture.py`); every reader goes through
`partial_pressures_pa`, so the storage could change underneath (since DB.71 it is the child tables). Mole fractions and the mean
molar mass follow from it. The rule that holds the stored values together:

- the partial pressures add up to `atmospheric_pressure_pa`;
- `atm_molar_density` is the mole-fraction mean of the gases' molar masses;
- `atm_density` is `P M / (R T)` and `scale_height_km` is `R T / (M g)` at the stored temperature and gravity;
- `surface_temperature_k` is the formula on the stored `bond_albedo` and `greenhouse_factor`;
- the ozone follows the oxygen and the star, and the aerosol depths follow their rules;
- every gas the text names is in the table and the mixture.

`tests/test_sweep_atmosphere.py` (TEST.209) checks these over planets of every class and zone round three host stars and the
bodies of generated systems, with the rules of `tests/audit/sweep_atmosphere.py` (`pressure-sum`, `mean-molar-mass`, `air-density`,
`scale-height`, `temperature`, `ozone`, `aerosol-airless`, `aerosol-cloud`, `aerosol-haze`, `aerosol-dust`, `text-gases`); its last test
breaks each rule in turn and checks that the sweep reports it.

**How the mixture is drawn.** `physics/atmosphere_recipes.py` draws the regimes (mantle redox, volatile delivery, outgassing,
retention, envelope, giant, brown dwarf, hot Jupiter, lava world, trace gases). `physics/mixture_draw.draw_mixture` turns them into
mole fractions by chemistry and not by a per-class template: the element budget, the temperature and pressure of the layer and
the oxygen fugacity go into the equilibrium solver (`physics/gas_equilibrium.py`, on the NASA Glenn thermodynamics of
`physics/species_thermo.py`). A class letter only gates which recipe rows apply. The family is chosen by what the recipe drew
(`mixture_rocky`, `mixture_giant`, `mixture_lava`, `mixture_escape`, `mixture_evolved`, and so on). The draws follow the recipe draw
in a fixed order; changing the order changes every body and is a generator epoch.

## 5. Ozone

`physics/ozone.py`. The column in Dobson units (1 DU = 2.687e16 molecules/cm2) is the fifth-degree polynomial of Kozakis et al.
2022 (A&A, Table 3) against the log10 of the oxygen mixing ratio, one for each of five hosts (G0V, Sun, K2V, K5V, M5V), valid for
0.01 to 150 percent of the present oxygen level and held at its ends outside it. The oxygen level at another pressure or gravity
is `(p_O2 / g)` over Earth's. Between two hosts the column is the linear mix in temperature; cooler than the M5V row takes it;
hotter than G0V scales the G0V curve to the F7V and F0V columns of [R3845, R3846] (a labelled stand-in). The draw multiplies
the median by a structure factor log-uniform on 0.85 to 1.71, which spans the published models. Below 0.01 percent of the present
level (abiotic oxygen) the column is 0.01 to 0.1 percent of the same host's full column, and 15 percent for an F host
([R1137], second-hand). No oxygen gives none. A loss event that kept a share of the air scales the column by that share
(`column_after_air_loss`).

## 6. Cloud, haze and dust

`physics/aerosols.py`, three things that scatter besides the gas:

- **Cloud.** The air is cooled along a dry adiabat to the skin temperature `T_s / 2^(1/4)`; a condensable gas forms a cloud where
  its partial pressure reaches saturation. The most condensate that can fall out gives the depth limit,
  `tau = 3 Q_ext m_c / (4 rho_p r_eff)` with `Q_ext = 2`; the droplet radius is log-uniform on 1 to 20 micron, the depth log-uniform
  from 0.1 to its limit, `cloud_cover` 0 without a condensate else uniform on 0.1 to 1. Clouds are grey, albedo 1, asymmetry 0.85.
- **Haze (stand-in).** Above a CH4/CO2 ratio of 0.1 methane photolysis makes a haze, supply limited by diffusion with `K_zz`
  log-uniform on 1e6 to 1e10 cm2/s and a yield log-uniform on 0.01 to 1; the depth is the particles' mass flux over their Stokes
  fall speed with the radius log-uniform on 0.05 to 0.3 micron and the Titan tholin density (1000 kg/m3). The extinction goes as
  lambda^-a with a from 0.97 to 2.34 and the absorption follows the Khare tholin index. The Titan laws are a test target, never a
  preset.
- **Dust (stand-in).** A dry solid surface lofts dust above the threshold friction velocity; `u*/u*_t` is log-uniform on 0.3 to 3, the
  depth log-uniform on 0.5 to 1.5, grey, radius 1.2 to 1.8 micron, asymmetry 0.673 +- 0.081, and the single-scattering albedo of
  the [R4743] spectra (digitised, a labelled stand-in).

The draws always use the same number of random numbers (four for the cloud, four for the haze, two for the dust, one for
`aerosol_draw`), whether or not the body has the aerosol, so one body's draws never shift another's.

## 7. Rayleigh scattering and the colour

`physics/rayleigh.py`: the optical depth of the scattering is `tau_R(lambda) = sum_i N_i sigma_i(lambda)`, with `N_i = p_i N_A / (m g)`
and `sigma = A lambda^-4 (1 + B / lambda^2)^2` fitted to the gas table's three values (450, 550, 650 nm). Dry Earth air gives 0.0970
at 550 nm against 0.09721 for [R0575: Table 4]. The cross sections are good to 2.5 to 4 percent for N2 and Ar, 7 for CO2 and 10
to 18 for CO, CH4 and N2O [R4138]. The colour model (`physics/atmosphere_colour.py`) is described in
[atmosphere-colour.md](atmosphere-colour.md); it reads the stored ozone column, the aerosol depths and the cloud cover, and a
stand-in ground albedo until the ground's reflectance is stored.

## 8. Limits

The infrared band and wing opacities are stand-ins until real line data are in (`physics/infrared_opacity.py`); the stored water
vapour is taken as it is, with no condensation above the surface; the closure's pressure exponent is 2 for every air. The haze yield
is a labelled placeholder (stand-in AT.n7) and the dust and haze spectra come from similar materials (stand-in AD2), both in
[stand-in-values.md](../design/research/stand-in-values.md); the ozone fit for a star hotter than G0V is labelled in the gap plan
([gap-filling-plan.md](../design/research/gap-filling-plan.md)).
