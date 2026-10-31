# Sky colour in planet and moon descriptions: analysis

**Status (2026-10-17):** the base of the atmosphere colour plan, [atmosphere-colour.md](atmosphere-colour.md). It was written as a read-only study of main at Release 8.0.1634 and is kept as filed; where Boss's later decisions differ (the albedo and greenhouse values are stored, the molar density is derived from the partial pressures, ozone, haze, dust and cloud are built in, gas giants are in, the colour is the absolute one) the plan note wins. "Sky colour" and "atmosphere colour" are the same thing (Boss, 2026-10-16, 11:24 PM EDT).


Read-only study of `dwhagar/planetGen` at main `1bb8755` (Release 8.0.1634). No code changed. Evidence tags: [C] read in the repo, [P] computed in a throwaway prototype (not in the repo), [S] source seen in a search result or fetched page, [R] recalled and not confirmed.

## Short answer

**Small to medium, and no migration is needed for a first version.** Almost everything a physical sky colour needs is already stored: the ten gas partial pressures, the surface pressure, gravity, the host star's temperature and luminosity, the orbital distance, the ozone-forming oxygen. The colour can be computed at description time as a pure function and shown as a new "Sky" row (hex swatch plus a colour word) in the planet and moon facts table. The two real gaps are **aerosols** (haze, dust, cloud: nothing is stored) and the **albedo** (drawn, used for temperature, thrown away). A first version without them is honest for clear, thin-to-moderate atmospheres and wrong for Mars-like dusty and Venus-like cloudy worlds; adding one stored column (albedo, plus an optional aerosol class) fixes most of that and does need a migration.

Prototype results (Rayleigh only, no ground, no aerosols) [P]: Earth air under the Sun gives a light blue (#9ac0ff) with an Earth Rayleigh optical depth of 0.098 at 550 nm, which matches the accepted 0.097 [R, Bucholtz 1995]. The same air under a 3,300 K red dwarf gives a cream-white sky in absolute terms, but the same blue once the eye is adapted to the star's light. A 1 MPa nitrogen world is white-cyan with an orange sun. A 9 MPa CO2 world is white with the sun disc extinguished. A thin clear CO2 atmosphere comes out blue, which real Mars is not (dust), so the aerosol gap matters.

## 1. What exists and what is missing

### Exists [C]

| Need | Where | Notes |
|---|---|---|
| Gas mix | `planets` and `moons` columns `p_o2_kpa p_co2_kpa p_co_kpa p_n2_kpa p_ar_kpa p_h2_kpa p_h2o_kpa p_ch4_kpa p_h2s_kpa p_so2_kpa` (`db/models.py:~740`, `physics/atmosphere.py`) | Partial pressures after the cold-side freeze-out. Set by GEN.85. |
| Total pressure | `atmospheric_pressure_pa` | Lowered by condensation. |
| Gravity | `gravity_g` | With P gives the air column: N = P / (m g). |
| Scale height | `scale_height_km` | Only set when there is air; computed from the airless equilibrium temperature, not the surface temperature (`physics/planets.py:~655-662`). Fine for a colour, but note it. |
| Surface temperature | `surface_temperature_k` | |
| Star | `stars.temperature_k`, `luminosity_w`, `radius_km`, `star_type`, `metallicity_fe_h`; `planet.distance` (AU) | A planet or moon reads its host through `planet.star` (`generation/planet.py` docstring, "star is the single source of truth"). A close binary is a `BinaryStarProxy` with `temperature` and `luminosity`. |
| Ozone proxy | `physics/radiation.py:216 ozone_ratio(p_o2_kpa)`, `uv_surface_index`, `ozone_loss_flag` | O3 is not stored as a gas, but its strength is already estimated from O2. |
| Surface type | `ocean_fraction`, `land_fraction`, `hydrosphere`, `ice_shell_km`, `life_chemical`, `reflection_spectrum_visible` (list of words) | Enough for a coarse ground and ocean tint. |
| Spin | `rotation_period_hours`, `axial_tilt_deg`, spin axis | Needed only for the sunset variant (a tidally locked world has no sunset). |
| Display path | `describe/build.py` (`_planet_facts`, `_atmosphere_composition`) builds a typed Document; `describe/render.py` makes text, markdown and wikitext; the web page and `GET /api/systems/<id>/description` consume the JSON (`docs/design/architecture/system-description.md`) | One place to add a row serves every output. |
| Related plan | `docs/design/orbits/sky-view.md` section 4 and stage 1 (`planetgen/physics/photometry.py`, not yet written) | Already specifies blackbody to sRGB through the Wyman, Sloan and Shirley (2013) CIE fit and an sRGB encoder. The sky colour reuses it; do not write it twice. |

### Missing

1. **Aerosols**: no haze, dust, cloud or sulfate-aerosol quantity anywhere (grep for haze, dust, cloud on bodies finds nothing). These decide Mars (dust), Venus (cloud deck), Titan (organic haze) and any wet world with cloud.
2. **Albedo**: `calculate_atmospheric_conditions` draws it from the class range (`tuning.PLANET_CLASSES[..]["albedo_range"]`, `physics/planets.py:657-658`) and keeps only the temperature it produced. It cannot be recovered afterwards because the greenhouse multiplier is also a random draw. The class range could stand in as a proxy (no migration), but a stored value is better and is also the best available cloudiness signal.
3. **Helium and the "balance" gas**: He (and ammonia, sodium) is named in the text but not stored (`atmosphere.py` `balance`). For H2/He giants He is about 14 percent by number; its Rayleigh cross section is about 70 times smaller than N2's, so ignoring it changes little. The mean molecular weight, however, should use the stored partials, not `atm_molar_density`, which is an independent random draw and does not match the partials.
4. **Ozone amount** and other absorbers' column densities. Only the O2 proxy exists.
5. **Sky-colour wording**: nothing says "sky" in a body description today. The only sky-like thing is the map glow, `glowColorForTemp` in `html/static/systemmap.js:129`, which is a three-colour temperature gesture, not physics. The sky colour could replace it later.

## 2. The simplest defensible method

All of this is a pure function `sky_colour(body)` in a new `planetgen/physics/sky.py`, called from `describe/build.py`. Nothing is fixed per class or copied from Solar System presets (project rule): every input is a stored physical quantity.

### 2.1 Light source

Planck spectrum at the host's `T_eff`, weighted by the CIE curves (reuse the sky-view stage-1 photometry code). Luminosity and distance only set the brightness (and whether the sun is bright enough for a "sky" at all); the colour depends on `T_eff`. Blackbody is adequate for F to K stars; for M dwarfs and cooler it ignores TiO and VO band structure, so state "blackbody approximation" in the doc and consider a model-atmosphere table later [R]. Binary: use the host star of the planet; a circumbinary planet uses the pair's luminosity-weighted mix.

### 2.2 Rayleigh scattering from composition and pressure

Optical depth per wavelength through the whole air column:

    tau_R(lambda) = N_col * sum_i f_i * sigma_i(lambda),   N_col = P / (m_bar * g)

with `f_i` the mole fractions from the stored partials, `m_bar = sum f_i m_i`, and the standard cross section

    sigma(lambda) = 24 pi^3 / (lambda^4 N_s^2) * ((n^2 - 1)/(n^2 + 2))^2 * F_K(lambda)

where `n` is the refractive index at standard density `N_s` and `F_K` the King correction for molecular anisotropy [S: formula quoted in papers citing Sneep and Ubachs 2005]. Per-gas `n(lambda)` and `F_K` come from the laboratory literature (below). Only ten gases, so a small constant table (Cauchy-type `n - 1 = A (1 + B / lambda^2)` per gas) is enough; no data files. The refractivities of N2, CO2, H2 and He are tabulated in Bétrémieux and Kaltenegger 2015 (arXiv:1507.02107, Table 1) [S]; the ratios to N2 that the prototype used came from those: CO2 about 2.5, H2 about 0.22, He about 0.014 (the sources check of 2026-10-17 corrects these to 2.5, 0.20 and 0.013; squared refractivity ratio times King factor) [P].

### 2.3 Sky radiance and the Sun's direct colour

- **Direct sun** (noon and any elevation): `S(lambda) = F_star(lambda) exp(-tau(lambda) * X)` where `X` is the airmass for the elevation. Noon is `X = 1`; sunset uses a Chapman-type function because the scale-height-to-radius ratio of the world sets the horizon airmass (`scale_height_km / radius_km` is already stored). Kasten and Young 1989 is the Earth form [R].
- **Blue sky, zenith or 90 degrees from the Sun**: single-scattering radiance with attenuation, `L(lambda) = F_star(lambda) * mu0/(mu0+mu) * (1 - exp(-tau (1/mu0 + 1/mu))) * (phase / 4 pi)`, with the Rayleigh phase function `3/4 (1 + cos^2 theta)` (Chandrasekhar 1960 [R]). This is the standard textbook limit and is accurate for `tau` up to about 0.3.
- **Thick air** (`tau_550` above about 0.5, roughly 5 bar of nitrogen and up): single scattering over-saturates the colour. Use a two-stream or adding-doubling pass (a short function, still no database) or fall back to the saturated limit: the sky tends to white or to the star's colour, the sun disc reddens, and at very large `tau` only the longest wavelengths reach the ground. The prototype shows this: 1 MPa N2 gives a white-cyan sky and an orange sun; 9 MPa CO2 gives a white sky and no sun disc [P].
- **Ground bounce** (surface albedo adds a floor of diffuse light, desaturating the blue): multiply by `1 / (1 - A_g * s)` with `s` the air's spherical albedo, Chandrasekhar's standard result [R]. Needs the ground tint (below).

### 2.4 Absorbers

- **Ozone** Chappuis band (500 to 700 nm) deepens the blue of Earth's sky. Estimate the column from `ozone_ratio(p_o2_kpa)` (already in `radiation.py`) times Earth's value, cross section from Burkholder and Talukdar 1994 or the Serdyuchenko 2014 compilation [R, not fetched]. A second-order effect; can come in phase 2.
- **Methane**: red and near-infrared bands tint ice-giant skies blue-green (Uranus, Neptune). Matters for gas-giant skies and for methane-rich small worlds; needs band cross sections [R]. Phase 2.
- **SO2, H2S, CO, H2O vapour**: weak in the visible; skip first.

### 2.5 Aerosols (Mie), the part with no data

Dust, haze and cloud scatter nearly grey (large particles) or absorb blue (iron-oxide dust, organic tholins), and they dominate Mars, Titan and Venus. A defensible minimum with no new data: derive an **aerosol optical depth class** from stored physics, not from the class letter.

- **Cloud deck**: from the stored albedo (a high Bond albedo with a thick atmosphere is cloud, as the repo's own class tuning notes) or from vapour supply: water vapour partial pressure against its saturation pressure at the surface temperature (the saturation function already exists as `vapour_pressure_kpa`). Near saturation gives cloud.
- **Dust**: scale with surface pressure, gravity and a dryness measure (low `ocean_fraction`, thin air raise dust lifting). No code or data for this today; it would be a new, openly flagged parametrisation.
- **Photochemical haze**: needs methane plus nitrogen plus UV (`uv_surface_index` exists). Titan-like. Flag as a later richer version.

Then add Mie extinction as a grey (slowly varying) optical depth plus an absorption slope for dusty or hazy air. Without the stored aerosol load the first version should say "clear-sky colour" explicitly.

### 2.6 Eyes and presentation

- **Adaptation is the main design choice.** Absolute colour under a red dwarf is cream or orange; a human adapted to that light would call it white or pale blue [P: 3,300 K gives #fffbea absolute versus #8cb6ff adapted]. Recommend showing the **adapted** colour (white-balanced to the star, a von Kries or Bradford step) as the swatch and naming the star's light in the sentence, with the absolute value available in the JSON. This is a presentation decision for Boss.
- **Ground and ocean tint** (the colour of the horizon band and the sea's reflected sky): ocean reflects the sky by Fresnel (about 2 percent at the zenith, rising towards the horizon); land colour from the same words `reflection_spectrum_visible` already carries for life worlds, or a neutral grey otherwise. Optional.
- **Words**: map the hue and lightness to a short phrase ("pale blue", "deep blue", "white", "pale yellow", "hazy orange", "dark") plus a hex.

## 3. Effort

| Version | What | Size |
|---|---|---|
| **1: clear-sky colour at noon** | `physics/sky.py` (Rayleigh table for ten gases, Planck to CIE to sRGB, single-scatter radiance, a thick-air limit, adaptation), one `fact("sky", "Sky", swatch, words)` row for bodies with air and a surface, a swatch inline item, unit tests. No DB change. | **Small to medium**, about 2 to 4 days of lane work including tests. |
| **2: aerosols, albedo and absorbers** | Store `bond_albedo` (planets and moons, one DOUBLE column, schema bump, backfill needs a regenerate because it cannot be recovered), aerosol class derived at generation time, Mie term, ozone and methane. | **Medium** for the code; the migration is a column add, and Boss's rule that the schema is not a constraint applies. |
| **3: sunset and twilight, swatch strip, gas-giant cloud-top sky** | Elevation-dependent colour, a day-to-night strip, optional map glow replacement. | **Medium**, naturally fits sky-view stage 6 (horizon view). |

Files that would change for version 1:

- new `src/planetgen/physics/sky.py` (and a shared `photometry.py` if sky-view stage 1 is built first or together);
- `src/planetgen/describe/build.py` (`_planet_facts`: one row, skipped for airless bodies, gas giants, and bodies without a star);
- `src/planetgen/describe/model.py` and `describe/quantity.py` (a new inline item `{"swatch": "#rrggbb", "name": "pale blue"}` that text, markdown and wikitext render as the word and the web renders as a chip; `walk_inline` and `check` need the new shape);
- the web HTML renderer for Documents (one new inline case) and `api.md` / `system-description.md` docs;
- tests: `src/tests/test_describe.py` (row present and absent cases, `check()`), new `test_sky.py` (Earth air gives blue, a red dwarf adapts to blue, tau scales with pressure, thick air saturates, no NaN for every class, deterministic), plus the existing describe parity tests.

No migration, no new table, no generator change (so no rescatter and no generator epoch bump) for version 1. The colour is deterministic from stored values, so it renders the same from freshly generated or database-loaded systems.

Test targets that are physical, not fixed content: Earth air column gives `tau_R(550 nm)` close to 0.097 [R, Bucholtz 1995 or Penndorf 1957]; the Rayleigh exponent is about -4.08 across the visible [R]; Earth's clear zenith sky is blue under a 5,772 K blackbody; a pure H2 atmosphere scatters about a fifth as much as N2 per molecule [S, arXiv:1507.02107 refractivities]; doubling pressure doubles `tau`.

## 4. Risks and a richer version

- **Aerosols decide the answer on many real worlds.** Version 1 on a thin CO2 world says blue (the prototype does), where Mars is butterscotch because of dust [S, search results: dust absorbs blue, so Rayleigh alone does not explain the colour] and its sunsets are blue from forward scattering by dust [S]. Say "clear-sky colour" in the label until aerosols exist, or hide the row below a pressure where dust would dominate.
- **Blackbody stars**: M dwarfs and hotter than about 10,000 K deviate; error is small for a swatch, larger for a deep red world.
- **Single scattering** is wrong for thick air; the plan above uses a limit or a two-stream step. Do not publish a precise hex for `tau_550` above 1 unless the multiple-scattering step is in.
- **Mean molecular weight**: use the stored partials, not `atm_molar_density` (they disagree because GEN.85 draws the mix independently of the older molar density). Worth a TODO for consistency whether or not sky colour ships.
- **Scale height input is the airless temperature** (`physics/planets.py:~655`). A real scale height changes only the sunset airmass, not the noon colour.
- **Perception model** (adaptation, rod/cone, gamut clipping): the sRGB gamut cannot show saturated violet or deep orange; normalisation hides the sky's absolute brightness. A "twilight to noon" strip needs an exposure model.
- **Gas giants and ice giants** have no surface; a "sky" is a pressure level (for example 1 bar). The prototype gave a pale blue at 1 bar for an H2/He column, but clouds and methane decide the real appearance. Suggest skipping them in version 1.
- **Moons of giants**: the primary fills the sky and reflects light; this is a different (and nice) feature, not part of sky colour.
- **Richer version**: noon colour, sunset colour and the sun-disc colour as three swatches; the horizon-versus-zenith gradient; a stored `bond_albedo` and aerosol load; ozone and methane bands; a full two-stream or Monte Carlo transfer; the sky-view horizon renderer (stage 6) using the same function per pixel. Compare with published full-atmosphere models (Bruneton and Neyret 2008, Hillaire 2020) [R] for Earth as a regression test.

## 5. Sources and how well each was checked

Checked in this session:
- Sneep and Ubachs, "Direct measurement of the Rayleigh scattering cross section in various gases", J. Quant. Spectrosc. Radiat. Transfer, 2005, DOI 10.1016/j.jqsrt.2004.07.025 [S]: cited in papers that quote the cross-section formula above. Volume and page (92, 293) not confirmed.
- Bucholtz, "Rayleigh-scattering calculations for the terrestrial atmosphere", Applied Optics 34(15), 2765-2773, 1995 [S]: citation confirmed from a reference list; content (tables and analytic fit from 0.2 to 1 micron) confirmed second-hand.
- Bétrémieux and Kaltenegger, "Refraction in planetary atmospheres: improved analytical expressions and comparison with a new ray tracing algorithm", arXiv:1507.02107 (ar5iv page fetched) [S]: refractivities at 1 micron N2 2.95e-4, CO2 4.43e-4, H2 1.37e-4, He 3.48e-5; H2, He.
- "The Habitability of Proxima Centauri b: II. Environmental States and Observational Discriminants", arXiv:1608.08620, a search result only (the full text was not retrievable; authors and journal not confirmed): says a planet's sky colour depends on atmosphere and starlight, might be pale blue or pale purple around Proxima, and that haze-covered Venus-like atmospheres look redder and lilac because Rayleigh scattering matters little in them. The sources check of 2026-10-17 read it: it is the colour of the planet's disk from space, a precedent for the star-plus-Rayleigh reasoning only, not for the colour of the sky from the ground.
- Mars: search results on dust absorbing blue and blue sunsets from forward scattering; an Applied Optics 2014 article (DOI 10.1364/AO.53.001808) is named in a result but its authors and title were not retrieved.
- Wyman, Sloan and Shirley (2013) CIE fit: already validated numerically in `docs/design/orbits/sky-view.md` [C].
- Environment note: arxiv.org is refused by the proxy; `ar5iv.labs.arxiv.org` served 1507.02107 but redirected to arxiv.org for 1608.08620, so that one was not read.

Recalled, not verified (check before citing in a design note): Peck and Huang 1977 (H2 refractivity, JOSA 67, 1550; a search did not return it), Penndorf 1957 (JOSA 47, 176), Chandrasekhar 1960 *Radiative Transfer*, Kasten and Young 1989 (airmass), Bruneton and Neyret 2008, Hillaire 2020, Burkholder and Talukdar 1994, Serdyuchenko et al. 2014 (ozone cross sections), Lemmon et al. 2004 Science 306, 1753 (Mars dust), Tomasko et al. 2008 and Karkoschka 1994 (Titan haze and ice giants).

## 6. Prototype (not committed)

The throwaway script is `/mnt/project-files/reports/sky-color-prototype-2026-10-17.py` (single scatter, Rayleigh only, relative cross sections from the table above, absolute 5.1e-27 cm2 for N2 near 532 nm scaled by lambda^-4.1). Numbers quoted above came from it; treat them as order-of-magnitude checks, not results. The 9 MPa sun-disc `#ff0000` printout is a clipping artefact of the single-scatter formula at tau of about 16 and is exactly the reason for the thick-air limit.
