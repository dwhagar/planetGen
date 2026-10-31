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
| surface temperature | 273 K in the paper; here the body's own when its climate is known, else 273 K | Table 2 | READ, body's own: GEN.735 |
| eta0 | uniform 9.3e7 to 1.9e11 Pa s (GEN.1185: the paper's 2.6e10 to 5.3e13 times 3.56e-3, so the viscosity at 1,600 K holds with the activation range below) | Section 2.2.3 retuned | stand-in IHA1 |
| Ea | uniform 260 to 390 kJ/mol (GEN.1185; the paper takes 200 to 300) | Korenaga and Karato 2008 Table 3, diffusion creep | READ, stand-in IHA2 |
| core radius | `R_dry CMF^0.5`; dry radius from `radius_fit_earth` or `small_body_radius_earth` | Zeng and Jacobsen 2017 via the paper; part B fits | READ |
| start | `Ts + E / C`, `E` = binding energy of the differentiated body (two uniform shells) | COMPUTED | COMPUTED |
| ceiling | the liquidus at the middle of the mantle, `a_l (p + p_l)^c_l` with 310.317, 29.102 GPa, 0.5337 (3,100 K for the Earth, 1,900 K for a Moon-sized body); a magma ocean radiates the rest away | Manske et al. 2025 (arXiv:2506.18718) Table 1, Eq. 7 | READ |

The two draws (`eta0`, `Ea`) are the paper's own uncertainty about other planets' mantles; the stream is two `uniform` calls. Everything else follows from the stored mass, radius, composition and the star's age.

Differences from the paper, all labelled here: one temperature instead of two; the stagnant lid's thickness is not followed (the lid is thin against the mantle, so its radiogenic heat counts with the mantle's); the integration starts at the body's formation from the computed start instead of priming for 5 Gy (the model is nearly memoryless, so the start matters only for a young body); the mantle is held at the 3500 K ceiling while the magma ocean and the first million years' 26Al would take it higher.

## 3. What it gives

| Case | Result |
|---|---|
| Earth twin, 4.568 Gy, paper's rheology range | 27 to 36 TW, Urey ratio 0.59 to 0.76, mantle 1400 to 2800 K, cooling 30 to 70 K per Gy |
| 0.1, 1 and 5 Earth masses at 4.5 Gy | Urey ratio 0.5 to 0.85 (the paper has about 0.65) |
| Earth twin at 0.1 Gy | 3 to 30 times the heat of 4.5 Gy, Urey ratio 0.1 to 0.25 (1.0 for the stiffest rheology, which stays near the ceiling) |
| Differentiation energy of the Earth | 1.9e31 J on top of 2.2e32 J of accretion |

The Earth's own measured heat flow is 42 to 47 TW with a Urey ratio of 0.47 to 0.52 (Nimmo et al. 2020; READ via part E 2.8e). This one box reaches neither with either lid (it holds the 21 TW of radiogenic heat and the cooling of one temperature); the totals are inside the geoneutrino range, 38.2 +13.6 -12.7 TW (Agostini et al. 2019).

### 3.1 The lid (GEN.721)

A body has a **tectonic regime**, `lid_regime`: `mobile` (plate tectonics) or `stagnant` (one rigid plate). It is drawn once, after the hydrosphere, at the chance the magnetic field's rule gave (`magnetism.mobile_lid_chance`: a liquid surface ocean and a mass of 0.5 to 5 Earth masses at 0.6, heavier at 0.3, any other body stagnant), and stored. Three things read it: the dynamo's lifetime (a stagnant lid's is 0.4 as long), the ground dose (a mobile lid reprocesses basalt into felsic crust) and the heat flux law here. Before GEN.721 the field and the ground dose each drew their own coin on the same chance, so one body could be mobile for its dynamo and stagnant for its dose.

The law. A stagnant lid loses heat across `dT = a_rh R Tm^2 / Ea`, the boundary layer's viscous temperature scale (section 1). A mobile lid carries the whole drop from the mantle to the surface, `dT = Tm - Ts`, with the same Rayleigh number and the viscosity at `Tm`, and a mobile exponent `beta = 1/3` (boundary-layer theory; Gurnis 1989 for a viscous lid with weak zones, READ in O'Neill 2020, arXiv:2006.05654 section 1; Moresi and Solomatov 1998 find 0.224 to 0.293 for a yielding lid and O'Neill's own runs 0.26 for an evolving mobile lid, 0.12 for a stagnant one). At one mantle temperature the mobile lid loses about ten times the heat, so the mantle cools to a lower temperature, the heat it was born with is spent sooner, and the radiogenic share is larger:

| Earth twin at 4.568 Gy, the paper's four rheology corners | stagnant lid | mobile lid |
|---|---|---|
| Heat out | 27 to 36 TW | 25 to 32 TW |
| Urey ratio | 0.59 to 0.76 | 0.66 to 0.83 |
| Mantle temperature (middle rheology) | about 1,950 K | about 1,430 K (the Earth's potential temperature is near 1,600 K) |

At the middle rheology and 0.1 Gy the mobile lid gives 1,007 TW against 701 TW; by 0.5 Gy the stagnant lid gives more. The ordering at 0.1 Gy depends on the rheology (at the low-viscosity corner the stagnant lid gives more).

The hydrosphere is made before the regime exists (the regime needs the hydrosphere), so it reads the heat of a stagnant lid; a surface ocean has no ice lid for that heat to thicken or thin, and the stored `internal_heat_flux_w_m2` is the regime's own, recomputed after the draw.

### 3.2 Small and large bodies, and a Urey ratio above 1 (GEN.763)

The Urey ratio is the radiogenic heat over the heat that leaves the surface (Stamenkovic et al. 2012, arXiv:1210.3838 section 3.1: "the ratio of internal heat production to surface heat flux"; they cite about 0.3 from geochemistry and close to one from the geophysical view for the Earth). Here it is an output, `radiogenic_w / total_w`. It passes 1 only while the mantle is warming (it makes more than it loses), and it is exactly 1 while the mantle sits at its ceiling.

Bulk silicate Earth composition, the middle rheology, a stagnant lid:

| Mass (Earth) | 0.1 Gy | 0.5 Gy | 1 Gy | 2 Gy | 4.5 Gy | 10 Gy |
|---|---|---|---|---|---|---|
| 0.0004 | 0.60 | 0.72 | 0.76 | 0.76 | 0.77 | 0.80 |
| 0.01 | 0.97 | 0.85 | 0.79 | 0.74 | 0.73 | 0.77 |
| 0.1 | 0.78 | 0.76 | 0.74 | 0.71 | 0.71 | 0.74 |
| 1 | 0.20 | 0.43 | 0.56 | 0.64 | 0.67 | 0.72 |
| 10 | 0.05 | 0.27 | 0.44 | 0.56 | 0.63 | 0.68 |

A small body made little heat to start with and has little left to give, so its radiogenic share is the larger; a large body keeps the heat of its assembly for longer, so it is the smaller. At 4.5 Gy, over both lids and the four rheology corners plus the middle, the ratio runs 0.49 (10 Earth masses, stagnant) to 1.0 (a stiff mantle held at its ceiling); at the middle rheology a mobile lid gives the larger ratio. The generated bodies of 24 sampled systems give medians of 0.78 (under 1e-4 Earth masses), 0.78, 0.76, 0.71, 0.67, 0.59 and 0.49 (2 to 5 Earth masses) in the mass bands 1e-4, 1e-3, 1e-2, 0.1, 0.5 and 2, from 0.28 to 1.0 in all but one case. The one above 1.0 is a physical warming: a moon of 0.0008 Earth masses around a 3.0 Gy star gives 1.25 with its mantle at 1,839 K and still warming, the radiogenic 24.6 GW against 19.7 GW out. A tidal body held at the heat-pipe temperature or the mantle ceiling sits at exactly 1.0.

The 3.0 and 4.8 that were seen for 0.0004 to 0.0007 Earth-mass moons were not the physics. The heat was stored when the body was made, at the star's provisional age; the star's age is raised afterwards (`Star.adjust_age_for_planets`), and the same body at its final age gives 0.71 to 0.76. The life pass in `StarSystem.__init__` now runs `interior_heat.update_heat` on every planet and moon once the age is final (no draws). What the other readers of the heat took before that (the hydrosphere's ice shell, the surface composition's volcanism gate) read the provisional age; they draw, so they cannot be redone on the run's stream. GEN.839 makes them again after the final-age heat (`planets.refresh_at_final_age`, called in the same pass): the hydrosphere with its ice shell, the ground with its volcanism gate, the ground dose when the ground changed and the habitability score, with their draws taken from the body's own stream (`draw.Stream`, seeded from what is final about the body), so no draw of the run's stream moves. The air is not made again: its degassing still reads the provisional age's heat (a deeper change, not part of this item).

### 3.3 The core's own temperature (GEN.838)

The mantle box above keeps one temperature for mantle and core, so the heat that crosses the core-mantle boundary, the inner core and the dynamo have no curve of their own. `physics/core_heat.py` adds a second state beside it, a core temperature at the top of the core, driven by the box's temperature history. The box is not changed by it: a body's stored heat flow and Urey ratio, and the mantle temperature, are what they were, and a body with no metal keeps no core (`Heat.core` is `None`).

What it says (a Stand-in throughout, the registry rows CC1 to CC11 in the gap plan section 6.136):

* The core is a uniform shell of the body's metal mass inside the core radius the box already uses, the mantle a uniform shell outside it. The pressures of the shells times 1.16 give the Earth's 136 GPa at the boundary and 337 GPa at the centre (364).
* The base of the mantle is the box's temperature plus 8 K per GPa of the pressure at the boundary. The core's heat crosses a boundary layer whose thickness comes from its critical Rayleigh number at the viscosity of the layer's mean temperature (the box's viscosity law plus a pressure term, the activation volume), so the heat is the drop over the layer. The core starts 1,000 K above the base of the mantle and cools as `C dT/dt = -Q`.
* The core's adiabat (6.74 K per GPa) meets the iron liquidus (anchored on the Earth's inner core boundary, 5,400 K at 329 GPa) where the inner core stops. The freezing iron adds its latent heat and gravitational energy per kelvin of cooling to the core's heat capacity, so the cooling slows.
* The dynamo runs while `(1 - T_b / T_c)(Q - Q_ad) + E_G dM/dt` is above zero, `Q_ad` the heat the core conducts down its own adiabat, and the liquid shell is at least a quarter of the core's radius.

For the Earth at 4.568 Gy (middle rheology) the core is at 4,080 K with an inner core of 0.40 of its radius (0.35 measured) and 5.7 TW crossing the boundary (5 to 15 TW), a dynamo of 2.4 TW. Venus keeps one, Mars has none (its dynamo stopped at 4.4 Gy in this model, the real one stopped within a Gy of formation: the model's weakest case), and a Mercury-like core freezes through (no light elements in the liquidus). A heavier planet's core stays liquid: at 2 to 10 Earth masses there is no inner core and the dynamo runs. `Core.dipole_scale` is the buoyancy-flux scaling of the moment over the Earth's, `(R / R_e)^2 (P / P_e)^(1/3) (rho / rho_e)^(1/6)`; `magnetism.generate_field` reads it, so the field comes from the core and the mantle's lid reaches it through the box's temperature (the lifetime draw of GEN.86 and the stagnant-lid factor of GEN.177 are gone).

The core is read where it is wanted (`interior_heat.core_of`), not stored: no column holds it. The field is drawn when the planet's spin is set, at the star's age then; a planet whose star is aged later (`adjust_age_for_planets`) keeps the field of that age.

## 4. Where it is used

`planets.internal_heat_flux_w_m2` (the hydrosphere's ice shell, `habitability_world`, the surface composition's volcanism gate) runs `interior_heat.update_heat` at the star's age at the time of the call and stores `urey_ratio` and `internal_heat_flux_w_m2` on the body; the mantle's two draws are stored once (`mantle_eta0_pa_s`, `mantle_activation_kj_mol`). `rogue_surface` runs the same model with a host-less content and the Earth-composition stand-in.

## 5. Stand-ins for the gap plan registry

* The separate core temperature (GEN.838, section 3.3): built as a second state driven by the box, not as the two-way coupled model of Thiriet et al. 2019 (the paper behind it, Physics of the Earth and Planetary Interiors 286, was not readable); the core's heat is not taken from the mantle's budget, so the mantle's history is Guimond et al. 2022's, with the core following it.
* The mantle's ceiling is the liquidus at the middle of the mantle (a uniform density and gravity); the solidus and a melt fraction are not modelled.
* The surface is the body's own surface temperature when its climate is known (`update_heat`), 273 K for a body without one (a mild effect, read only by a mobile lid's drop and by the birth temperature); a mantle colder than `Ea / R T` = 600 does not convect and loses no heat (no conduction), so a cold small body warms with its radiogenic heat until convection starts.
* The mobile lid's exponent (1/3) and critical Rayleigh number (450, the stagnant lid's) are unsourced for this model; replace them with a calibrated mobile-lid law (Korenaga, Driscoll and Bercovici) when one is read.
* The regime is a coin on a mass and ocean rule (GEN.177), not a function of the thermal history; whether a body has plate tectonics is unsettled (the survey H2).
