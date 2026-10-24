# Interior heat: radiogenic plus primordial, the Urey ratio an output (GEN.699)

Status: built 2026-10-24 by lane H on top of GEN.459 (`physics/radiogenic.py`). Code: `src/planetgen/physics/interior_heat.py`; tests: `src/tests/test_interior_heat.py`; two reference checks in `physics/mathcheck.py`. Boss (1:28 AM EDT 2026-10-24, yes to the coordinator's question): model the heat a body is born with separately, so the total is radiogenic plus primordial and no Urey ratio is drawn.

## 1. What changed

GEN.459 drew a Urey ratio (0.1 to 0.9) and divided the radiogenic heat by it. That hid the formation heat and made a young body give the same share of radiogenic heat as an old one. Now a body has one thermal history. The radiogenic term is the stored content of GEN.459 (K, Th, U and 26Al at formation, decayed to the age); the primordial term is the cooling of the interior, and the Urey ratio is the first over the sum.

## 2. The model

One-box parameterised convection for a stagnant lid: Thiriet et al. 2019 as used by Guimond, Rudge and Shorttle 2022 (arXiv:2201.05636 section 2.2 and Table 2; READ). Mantle and core share one temperature `Tm` (that paper's runs have the core following the mantle):

`(Mm cm + Mc cc) dTm/dt = -qu Au + Mm h(t)`

`qu = km dT / delta`, `delta = d (Ra_crit / Ra)^(1/3)`, `Ra = alpha rho g dT d^3 / (kappa eta(Tm))`, `dT = a_rh R Tm^2 / Ea`, `eta = eta0 exp(Ea / R Tm)`.

| Quantity | Value | Source | Label |
|---|---|---|---|
| rho, cm, cc, km, alpha, kappa | 3500 kg/m3, 1142 and 840 J/kg/K, 4 W/m/K, 2.5e-5 /K, 1e-6 m2/s | Guimond 2022 Table 2 (from Thiriet 2019) | READ |
| Ra_crit, a_rh, beta | 450, 2.44, 1/3 | Table 2 (Thiriet 2019; Solomatov 1995) | READ |
| surface temperature | 273 K | Table 2 | READ |
| eta0 | uniform 2.6e10 to 5.3e13 Pa s | Section 2.2.3 (wet olivine less a decade to dry enstatite) | READ |
| Ea | uniform 200 to 300 kJ/mol | Section 2.2.3 (Karato and Wu 1993; Zhang et al. 2017) | READ |
| core radius | `R_dry CMF^0.5`; dry radius from `radius_fit_earth` or `small_body_radius_earth` | Zeng and Jacobsen 2017 via the paper; part B fits | READ |
| start | `Ts + E / C`, `E` = binding energy of the differentiated body (two uniform shells) | COMPUTED | COMPUTED |
| ceiling | 3500 K (a magma ocean radiates the rest away) | the liquidus of a deep magma ocean | STAND-IN |

The two draws (`eta0`, `Ea`) are the paper's own uncertainty about other planets' mantles; the stream is two `uniform` calls. Everything else follows from the stored mass, radius, composition and the star's age.

Differences from the paper, all labelled here: one temperature instead of two; the stagnant lid's thickness is not followed (the lid is thin against the mantle, so its radiogenic heat counts with the mantle's); the integration starts at the body's formation from the computed start instead of priming for 5 Gy (the model is nearly memoryless, so the start matters only for a young body); the mantle is held at the 3500 K ceiling while the magma ocean and the first million years' 26Al would take it higher.

## 3. What it gives

| Case | Result |
|---|---|
| Earth twin, 4.568 Gy, paper's rheology range | 27 to 37 TW, Urey ratio 0.57 to 0.76, mantle 1400 to 2800 K, cooling 30 to 70 K per Gy |
| 0.1, 1 and 5 Earth masses at 4.5 Gy | Urey ratio 0.5 to 0.85 (the paper has about 0.65) |
| Earth twin at 0.1 Gy | 4 to 33 times the heat of 4.5 Gy, Urey ratio 0.1 to 0.2 (0.6 for the stiffest rheology, which stays near the ceiling) |
| Differentiation energy of the Earth | 1.9e31 J on top of 2.2e32 J of accretion |

The Earth's own measured heat flow is 42 to 47 TW with a Urey ratio of 0.47 to 0.52 (Nimmo et al. 2020; READ via part E 2.8e). A mobile lid (plate tectonics) loses heat more freely than a stagnant one, so an Earth twin on this model is under it; it is inside the geoneutrino total, 38.2 +13.6 -12.7 TW (Agostini et al. 2019). Whether a body has a mobile lid is drawn after the hydrosphere (`magnetism.mobile_lid_chance`), so a tectonic-mode switch of the heat loss is a separate item.

## 4. Where it is used

`planets.internal_heat_flux_w_m2` (the hydrosphere's ice shell, `habitability_world`, the surface composition's volcanism gate) runs `interior_heat.update_heat` at the star's age at the time of the call and stores `urey_ratio` and `internal_heat_flux_w_m2` on the body; the mantle's two draws are stored once (`mantle_eta0_pa_s`, `mantle_activation_kj_mol`). `rogue_surface` runs the same model with a host-less content and the Earth-composition stand-in.

## 5. Stand-ins for the gap plan registry

* Ceiling temperature 3500 K (`tuning.INTERIOR_INITIAL_TEMPERATURE_CAP_K`): the liquidus of a deep magma ocean; replace with a pressure-dependent liquidus when the mantle's melting is modelled.
* The surface held at 273 K for the boundary layer, as the paper does; replace with the body's surface temperature once heat flow and climate are coupled (a mild effect).
* Stagnant lid only; a mobile-lid cooling law waits on the regime.
