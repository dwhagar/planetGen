# Atmosphere retention, hot- and cold-zone classes and life-stage timing

When a planet can hold an atmosphere, what the hot-zone and cold-zone
counterparts of Class S and V look like, how the planet classes should be
re-cut around that, and how the highest life stage should follow the time a
planet has been habitable. This is the research behind `docs/TODO.md` items
GEN.91 and GEN.92 and the class work they gate (GEN.28, GEN.27, GEN.29,
GEN.33, GEN.38). It carries the still-valid formulas of Boss's uploaded
documents "Mathematical and Algorithmic Implementation of the Planetary
Habitability Index.md" (Math), "Atmospheric Toxicity.md" (Atm) and "Planetary
Habitability and Speculative Xenobiology.md" (Master), and names every place
they are corrected (section 8). The score structure and the equipment tiers
are in [habitability-index.md](habitability-index.md); stellar activity,
magnetism, radiation dose and oceans are in
[activity-magnetism-radiation-hydrosphere.md](activity-magnetism-radiation-hydrosphere.md).

Informs: GEN.91, GEN.92, GEN.85, GEN.90, GEN.28, GEN.27, GEN.29, GEN.33, GEN.38
Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [S] seen in a search result (URL under Sources), [C]
computed (formula and inputs shown, reproducible), [R] recalled from memory
and unconfirmed, [D] a design default with no outside source. The research
environment could only read search-result text, not the papers, so every [R]
is listed under Evidence notes for checking when paper access is allowed.

## Decisions already taken

- Boss (2026-10-03 05:38Z), GEN.91: "We need a planet class for like Class S
  and V but in the other zones, so that will require research on when an
  atmosphere is possible and when it isn't."
- Boss (2026-10-03 05:38Z), GEN.90: "Refactor all planetary classes to better
  align with our habitability index system and the research for the
  different kinds of life, etc."
- Boss (2026-10-07 17:11Z, confirmed 2026-10-09 07:48Z): the score structure
  is PHI-4's domains and tiers for display with the Xenobiology doc's three
  scores behind them ([habitability-index.md](habitability-index.md)); GEN.29
  stays in phase 2 with the class refactor.
- GEN.33's rule: one planet class per PR, each with its tests.

## 1. Summary

1. **An atmosphere is possible when three gates pass**, each a pure function of
   mass, radius, insolation, star age and surface temperature (code in 3.6).
   - Cosmic shoreline: cumulative-XUV-weighted insolation below
     `S_cs = 3.85e-4 * v_esc^4` (S in Earth units, v_esc in km/s). Ratio
     `r = I_xuv / S_cs`: under 0.7 dense, 0.7 to 1 moderate, 1 to 2.5 tenuous
     (Mars), 2.5 to 30 trace or outgassed only, above 30 none.
   - Species: a gas is held thermally if `lambda = m v_esc^2 / (2 k T_exo)` is
     at least 30.
   - Cold side: a gas below its vapour-pressure temperature is ice, not air.
   - H/He envelopes and oceans also need the energy-limited XUV check.
2. **Answer to Boss's GEN.91 question.** The "S and V in the other zones"
   classes are mostly Class S itself with air, plus a hot desiccated class
   (Venus type, a re-scoped N) and a cold ice-shell class (an extended X, or a
   Hycean flag on R). Class S has `atmosphere = None` in every zone, but a 2 to
   10 Earth-mass rock in the ecosphere or cold zone sits far inside the
   shoreline (r about 0.02) and is airless only above T_eq of about 650 to
   1,200 K [C].
3. **No new class letters are needed** (5.2): A to Z are full after R, U, W, X,
   Y and Z, and Q may exist in saved galaxies.
4. **Life stage follows time habitable and the score gates** instead of a
   per-class cap (section 7).

## 2. What the code does today (checked 2026-10-09)
- **Classes.** `tuning.PLANET_CLASSES` has zone flags `h/e/c/r`,
  `radius_range`, optional `density_range`, an `atmosphere` text label and
  per-class albedo and greenhouse ranges. No retention test exists.
- **Pressure** is `atm_density * g * H` times `_atmosphere_retention_factor(g)
  = g^1` (`physics/planets.py`, line 580), "a starting point, not derived".
  Temperature, XUV and age play no part.
- **Habitable zone.** `orbits.calculate_habitable_zone` uses fixed
  `S_inner = 1.1` and `S_outer = 0.53` for every star: for the Sun 0.95 and
  1.37 AU against Kopparapu's 0.95 (runaway) and 1.68 AU (maximum greenhouse)
  [C].
- **Stellar evolution.** `evolve_star` holds a main-sequence star at constant
  luminosity: no pre-main-sequence phase and no brightening.
- **Life.** `HABITABLE_PLANET_CLASSES = E F G H K L M O P V`;
  `PLANET_CLASS_MAX_LIFE_STAGE` caps E (abiogenesis), F and G (photosynthesis)
  and L (multicellularity). `evolution.get_evolutionary_timeline` picks the
  latest milestone whose absolute age is below the star's age. Nothing reads
  the planet's score, atmosphere, radiation or time since habitable.
- **A finding.** Class N (and Q) carries a `life_chemical` but is not in the
  habitable list and has no stage cap. `life.py` (line 198) stores an uncapped
  `evolutionary_data` timeline for any ecosphere planet with a chemical;
  `Planet.to_paragraph_list` displays it only for the habitable list. A Venus
  analog (737 K, 9.2 MPa) therefore stores a life chemical and a hidden
  timeline, and its page still says "conditions are suitable for the
  development of life".
- **Tidal locking** exists only for moons; a planet's spin is a uniform 10 to
  1,400 h draw (GEN.104 adds spin). `planet_class` is `String(16)` on planets
  and moons, `String(4)` on `rogue_planets`.

Generator medians around a G2V star (40 planets per class and zone) [C]:

| Class, zone | Median R (km) | M (Earth) | v_esc (km/s) | S (Earth) | T (K) | P |
|---|---|---|---|---|---|---|
| N ecosphere | 6,072 | 0.77 | 10.0 | 1.05 | 743 | 9.4 MPa |
| S hot / e / c | 8,810 to 8,880 | 3.7 to 3.9 | 18.4 | 4.4 / 0.74 / 0.13 | 376 / 240 / 157 | 0 |
| V ecosphere | 11,245 | 4.9 | 18.8 | 0.88 | 382 | 268 kPa |
| I ecosphere / T cold | 26,200 to 27,700 | 20 to 22 | 23 to 25 | 0.74 / 0.13 | 242 / 159 | about 0.2 MPa |

I (10 to 60 Earth masses) and T (6 to 40) overlap almost completely.

## 3. Atmosphere retention

### 3.1 Gate 1: the cosmic shoreline
Zahnle and Catling (2017, arXiv:1702.03386): solar-system bodies with and
without atmospheres separate along cumulative insolation proportional to
`v_esc^4`, and exoplanets crowd along the extrapolation [S]. A pure
energy-limited argument gives an exponent of about 3, so the slope is empirical
[S]. Meni-Gallardo and Pallé (2025) keep slope 4 with a higher zero point [S,
snippets only]. The constant here is fitted, not copied.

Fit [C], bodies as (v_esc km/s, S in Earth units, state): Mercury (4.25, 6.67,
none), Venus (10.36, 1.91, dense), Earth (11.19, 1.0, dense), Moon (2.38, 1.0,
none), Mars (5.03, 0.431, 6 mbar), Ganymede (2.74, 0.037, none), Callisto
(2.44, 0.037, trace), Io (2.56, 0.037, trace), Europa (2.03, 0.037, none),
Titan (2.64, 0.011, 1.5 bar), Triton (1.45, 0.0011, 14 micro-bar), Pluto (1.21,
0.00064, 10 micro-bar), Ceres (0.51, 0.13, none). For
`S_cs = K v^4`, Titan forces `K > 2.26e-4`, Pluto and Triton raise the lower
bound to 3.0e-4, and Ganymede forces `K < 6.56e-4`. The window is 3.0e-4 to
6.6e-4; the geometric centre of the wider window, 3.85e-4, is used. Ratios `r`
at K = 3.85e-4: Earth 0.17, Venus 0.43, Titan 0.59, Triton
0.65, Pluto 0.78, Ganymede 1.7, Mars 1.8, Io 2.2, Callisto 2.7, Europa 5.7,
Mercury 53, Moon 82, Ceres 4,990.

| v_esc (km/s) | S_cs (Earth = 1) | T_eq at S_cs (K, A = 0.3) |
|---|---|---|
| 2 | 0.006 | 71 |
| 4 | 0.099 | 143 |
| 6 | 0.50 | 214 |
| 10 | 3.9 | 357 |
| 15 | 19.5 | 536 |
| 20 | 62 | 714 |
| 30 | 312 | 1,071 |

**M dwarfs need cumulative XUV, not bolometric insolation**: they stay
saturated for billions of years. Replace S by `I_xuv = S * (cumulative XUV of
the star to its age / cumulative XUV of the Sun to 4.5 Gyr)`, the "XUV exposure
index"
([activity-magnetism-radiation-hydrosphere.md](activity-magnetism-radiation-hydrosphere.md)
2.3). With a stand-in XUV history (`t_sat` 2 Gyr, decay index 1.23) TRAPPIST-1 b and
c come out at r = 5.2 and 3.0 (JWST: bare rock or under 0.1 bar [S]), e at 2.2,
f at 0.71 and g at 0.35 (could hold air); LHS 3844 b at 36 (none; bare rock
[R]); 55 Cnc e at 22 (JWST finds a likely CO2/CO atmosphere over a magma ocean
[S, Hu et al. 2024], so the 2.5 to 30 band is real).

The shoreline caps pressure, it does not set it: Titan, Triton and Pluto are
retention-allowed but their pressures come from supply and vapour pressure.
The generator takes `P = min(P_supply, P_retention, P_vapour)`.

### 3.2 Gate 2: which gases are held
Jeans parameter `lambda = G M m / (k T_exo r_exo) = m v_esc^2 / (2 k T_exo)`
(`r_exo` set to R). Kept if `lambda >= 30`, lost within about 1 Gyr if under 15
[R: Catling and Kasting 2017]. A lifetime integral for a 1-bar pure-gas
inventory gives Gyr or longer for `lambda_exo` of about 10 or more [C], so 30 is
conservative.

Exobase temperature is the hard input. A one-line model that fits Earth
(1,004 K), Venus (389 K model; about 275 K observed day side), Mars (268 K),
Titan (162 K) and Pluto (60 K model, about 70 K observed) is [C, a calibration,
not physics]:

```
T_exo = 1.3 * T_eq                                        for CO2-dominated air
T_exo = T_eq + min(1500, 750 * sqrt(F_xuv_rel))   K       otherwise
```

`F_xuv_rel` is the star's present XUV flux at the planet in units of Earth's
today (4.6 erg/cm2/s; the activity document's table 2.2).

Minimum retained molecular weight `mu_min = 2 k T_exo * 30 / (amu v_esc^2)`
[C]; a lighter gas is lost over Gyr. Weights: H2 2.0, He 4.0, CH4 16, H2O 18,
N2 and CO 28, O2 32, Ar 40, CO2 44, SO2 64.

| v_esc (km/s) \ T_exo (K) | 200 | 500 | 1,000 | 2,000 |
|---|---|---|---|---|
| 3 | 11.1 | 27.7 | 55 | 111 |
| 5 | 4.0 | 10.0 | 20.0 | 39.9 |
| 10 | 1.0 | 2.5 | 5.0 | 10.0 |
| 15 | 0.44 | 1.1 | 2.2 | 4.4 |
| 20 | 0.25 | 0.62 | 1.25 | 2.5 |
| 30 | 0.11 | 0.28 | 0.55 | 1.1 |

Earth (11.2 km/s, 1,000 K) holds everything above about 4: H2 and He are lost,
as observed. Planets of 6 to 10 Earth masses keep H2 against Jeans escape below
T_exo of about 500 K; XUV loss then decides (3.4). Jeans alone would keep H2 on
a Venus-like planet (10.4 km/s, 390 K), which is wrong because hydrogen leaves
by hydrodynamic and non-thermal routes: the species gate is necessary, not
sufficient.

**Water needs the H atom.** Water is lost by photolysis plus escape of atomic
H. Mars lost its water because H escapes; Earth keeps its ocean because the
cold trap leaves the stratosphere dry. The water check is the cold-trap
condition (4.1) followed by the energy-limited H clock (3.4).

**Non-thermal escape** (photochemical, pick-up, sputtering) explains Mars and
the Moon [R]. A global field is second-order for a thick atmosphere, so the
shoreline has no field term.

### 3.3 Gate 3: the cold side

Minimum temperature (K) for a gas to be vapour at partial pressure P, from
Clausius-Clapeyron anchored at the 1-bar boiling or sublimation point [C].
Checks: CO2 at 0.6 kPa gives 148 K (Mars polar frost point); N2 at 1 Pa gives
33 K against about 37 K at Pluto.

| gas \ P (kPa) | 100 | 10 | 1 | 0.1 | 0.01 | 0.001 |
|---|---|---|---|---|---|---|
| H2 | 20 | 14 | 11 | 9 | 7 | 6 |
| N2 | 77 | 61 | 50 | 43 | 37 | 33 |
| CO | 81 | 65 | 54 | 46 | 40 | 36 |
| Ar | 87 | 69 | 57 | 49 | 43 | 38 |
| O2 | 90 | 72 | 60 | 51 | 45 | 40 |
| CH4 | 112 | 88 | 73 | 63 | 55 | 48 |
| CO2 | 195 | 169 | 150 | 135 | 122 | 112 |
| NH3 | 240 | 200 | 172 | 151 | 134 | 121 |
| SO2 | 263 | 219 | 187 | 164 | 145 | 131 |
| H2O | 373 | 317 | 276 | 253 | 231 | 212 |

A planet with T_surface below the entry at its target partial pressure holds
that gas as ice (Titan at 94 K holds 1.5 bar of N2 and condensing CH4; Pluto at
40 K only about 1 Pa of N2). A CO2 or N2 atmosphere on a locked planet
collapses onto the night side unless it is about 0.1 bar or more [R: Joshi et
al. 1997]. The rogue-planet model already uses the N2 version
(`rogue_surface.nitrogen_vapor_pressure_pa`).

### 3.4 Energy-limited XUV loss: envelopes and oceans
`Mdot = eps * pi * R^3 * F_xuv / (G M)` with eps = 0.1, `K_tide` = 1 and
`R_XUV = R`, a lower bound on loss (Math's section 2.1 formula carries `R_p
R_XUV^2` and `K_tide`). Integrated over the star's life it is a fraction of
planet mass `f_lost = eps pi R^3 E_cum / (G M^2)`, with `E_cum` the cumulative
XUV energy per area. It is an upper bound (every absorbed photon lifts gas),
and for water it also needs a wet stratosphere.

Insolation at which the cumulative loss equals a fraction of a rocky core's
mass, Sun-like star (`E_cum = 3.1e16 J/m2` at S = 1 over 4.5 Gyr, from the
activity document's XUV model [C]; core radius `M^0.27`):

| M_core (Earth) | S for 1 percent | S for 3 percent | S for 10 percent |
|---|---|---|---|
| 1 | 9 | 28 | 93 |
| 2 | 21 | 64 | 213 |
| 5 | 63 | 190 | 633 |
| 10 | 144 | 433 | 1,443 |
| 20 | 329 | 988 | 3,293 |

A saturated fraction of 1e-3 (the X-ray band alone) would give S values about
twice as large. Radius-valley planets orbit Sun-like stars at S of 100 to 300 [R], where
envelopes of 1 to 3 percent of mass are stripped either way. Around an M dwarf
(`t_sat` 2 Gyr) the same loss happens at about a tenth of the insolation.

Ocean loss clock [C]: the H in one Earth ocean (1.55e20 kg) is removed in 14 to
25 Myr at S = 1.1 and 5 to 9 Myr at S = 3 during the saturated phase, if loss
is energy-limited and the stratosphere wet; 100 Myr of runaway greenhouse costs
4 to 17 oceans' worth of H. Luger and Barnes found several oceans lost and
hundreds to thousands of bar of abiotic O2 for habitable-zone planets of M
dwarfs older than about 1 Gyr [S]. Desiccation rule: `water_lost` = sum over
time-in-runaway of `Mdot_H`; desiccated if it exceeds the ocean inventory
(default one ocean; a log-normal around 1 to a few oceans is better).

### 3.5 Sub-Neptunes and the radius valley

- The radius valley (Fulton et al. 2017) is a factor-2 deficit near 1.5 to 2.0
  Earth radii [S]. Owen and Wu (2017) explain it by photoevaporation: bare
  cores near 1.3 Re, envelope planets near 2.6 Re, sculpted in the first
  100 Myr [S]; its period slope is negative, `R_valley(P) = 1.8 Re (P/10 d)^-0.09`
  [R]. Core-powered mass loss gives the same observable [S].
- Chen and Kipping (2017): Terran-to-Neptunian break at 2.0 (+0.7, -0.6) Earth
  masses [S]; median radii 1 Me 1.01 Re, 3 Me 1.54, 5 Me 2.07, 10 Me 3.11,
  20 Me 4.67 (slope 0.279 below 2.06 Me, 0.586 above) [C].
- Earth-like composition `R = M^0.27` [R: Zeng et al.] gives 6 Me = 1.62 Re.
  Classify on bulk density, not radius: rocky if `R <= 1.12 M^0.27` [D],
  volatile-rich (class R) otherwise. The repo's `GIANT_NEPTUNIAN_MASS_RADIUS`
  (`physics/constants.py`, line 179) gives 6 Me 2.26 Re and 20 Me 4.4 Re.

### 3.6 A generator-evaluable function

Tested on the bodies of 3.1 (Earth r 0.17 dense, Venus 0.43, Mars 1.77
tenuous, Mercury 53.7 none, Moon 81.6 none, Titan 0.59) [C]. `xuv_index` and
`xuv_rel_now` come from the activity module (pass S for both for
solar-system bodies).

```python
import math
AMU, KB, G, ME, RE = 1.66054e-27, 1.380649e-23, 6.674e-11, 5.972e24, 6.371e6
K_SHORE, LAMBDA_KEEP = 3.85e-4, 30.0            # K window 3.0e-4 .. 6.6e-4
MU = {'H': 1.008, 'H2': 2.016, 'He': 4.003, 'CH4': 16.04, 'H2O': 18.015, 'N2': 28.01,
      'O2': 32.0, 'Ar': 39.95, 'CO2': 44.01, 'SO2': 64.07}

def v_esc_kms(m_e, r_e): return 11.186 * math.sqrt(m_e / r_e)

def retention(m_e, r_e, s, xuv_index, xuv_rel_now, albedo=0.3, co2_dominated=False, k=K_SHORE):
    v = v_esc_kms(m_e, r_e)
    r = xuv_index / (k * v ** 4)                          # shoreline ratio
    teq = 278.6 * (s * (1 - albedo)) ** 0.25
    texo = 1.3 * teq if co2_dominated else teq + min(1500.0, 750.0 * math.sqrt(xuv_rel_now))
    kept = [g for g, mu in MU.items() if mu * AMU * (v * 1e3) ** 2 / (2 * KB * texo) >= LAMBDA_KEEP]
    verdict = ('dense' if r < 0.7 else 'moderate' if r < 1 else 'tenuous' if r < 2.5
               else 'trace' if r < 30 else 'none')
    return v, r, teq, texo, verdict, kept

def envelope_loss_fraction(m_e, r_e, cum_xuv_j_m2, eps=0.1):
    # cum_xuv_j_m2: cumulative XUV energy per area at the planet to the planet's age
    return eps * math.pi * (r_e * RE) ** 3 * cum_xuv_j_m2 / (G * (m_e * ME) ** 2)
```

Draw `k` log-uniformly from 3.0e-4 to 6.6e-4 per planet (or per system) to
reproduce the scatter between Titan and Ganymede. `P_retention` is the only
new piece of the final pressure `min(P_supply, P_vapour(T_surf, gas),
P_retention(r))`, and can be a smooth cap: `r < 0.7` none, 0.7 to 1: 0.1 bar,
1 to 2.5: 10 mbar, 2.5 to 30: 10 micro-bar, above: 0.

## 4. Habitable-zone limits, hot-zone outcomes and Hycean worlds

### 4.1 Runaway greenhouse and zone edges

Kopparapu et al. (2013, ApJ 765:131; 2014, ApJL 787:L29):
`S_eff(Teff) = S0 + a T* + b T*^2 + c T*^3 + d T*^4`, `T* = Teff - 5780`,
valid 2,600 to 7,200 K. Confirmed [S]: for the Sun recent Venus 1.78, runaway 1.04 (2013) or 1.107 (2014,
1 Earth mass), moist 1.01, maximum greenhouse 0.35, early Mars 0.32; 2014
runaway a = 1.332e-4, b = 1.58e-8; maximum greenhouse S0 0.356, a = 6.171e-5,
b = 1.698e-9; mass dependence only for the inner edge. The c, d coefficients
are [R], but reproduce the quoted solar distances (0.75, 0.95, 0.993, 1.676,
1.774 AU) [C]. The 5 Earth-mass runaway set (S0 1.188) and 0.1 Earth-mass set
(S0 0.99) are [R] and not tabulated; take them from Kopparapu 2014.

| Boundary | S0 | a | b | c | d |
|---|---|---|---|---|---|
| Recent Venus | 1.776 | 2.136e-4 | 2.533e-8 | -1.332e-11 | -3.097e-15 [R] |
| Runaway, 1 Earth mass | 1.107 | 1.332e-4 [S] | 1.58e-8 [S] | -8.308e-12 [R] | -1.931e-15 [R] |
| Moist (2013) | 1.014 | 8.1774e-5 | 1.7063e-9 | -4.3241e-12 | -6.6462e-16 [R] |
| Maximum greenhouse | 0.356 | 6.171e-5 [S] | 1.698e-9 [S] | -3.198e-12 [R] | -5.575e-16 [R] |
| Early Mars | 0.3179 | 5.4513e-5 | 1.5313e-9 | -2.7786e-12 | -4.8997e-16 [R] |

Runaway `S_eff` for 1 Earth mass: 2,800 K 0.918; 4,000 K 0.947; 5,780 K 1.107;
7,200 K 1.296 [C]. The repo's fixed 1.1 and 0.53 put the inner edge of a
3,200 K M dwarf 8 percent too far out in AU and its outer edge (0.53 against
0.238) 33 percent too near [C]; for a Sun-like star the inner edge matches and
the outer is 1.37 AU against 1.68. Recommended: a `hz_limits(teff, lum,
planet_mass)` function with this table, as its own item (it changes zone
assignment for every star).

Post-runaway history matters: a planet above `S_runaway` early (pre-main-
sequence M dwarf) and desiccated stays dry after the star settles, and the
repo has no pre-MS luminosity. Minimum model: for stars under about 0.5 Msun
give planets in the present zone an extra runaway time of 0.1 to 1 Gyr [R;
Luger and Barnes: "several hundred Myr" [S]]; for Sun-like stars use Gough's
`L(t)/L_now = 1/(1 + 0.4 (1 - t/t_now))` rescaled by `t/t_MS` [R].

### 4.2 Hot-zone outcomes
- **Venus type (desiccated).** `S_runaway(Teff, M_p) < S < S_cs(v_esc)` plus
  enough XUV loss (3.4): CO2 of tens of bar (Venus 92), N2 3 to 4 bar,
  SO2/H2SO4 clouds, T_surface 450 to 900 K, no standing water, no
  carbonate-silicate thermostat. Beyond `S_cs` the planet is bare lava rock
  (Class S hot). This is the default hot subtype. Atm and Master describe a
  1-bar "dune world" with masks and evaporative cooling (the Abe et al. 2011
  dry-land planet, stable in a narrow band); keep it as a rare variant (about
  10 percent).
- **Lava / magma-ocean world.** T_eq above about 1,500 K (55 Cnc e about
  2,000 K [S]): SiO, Na, K, O vapour if r is above 30, magma-ocean CO2/CO if
  2.5 to 30.
- **Steam worlds** (water over 10 percent, very close in) are not stable for
  Gyr; leave them out. **Stripped cores:** hot sub-Neptunes after
  photoevaporation are Class S; R and S are the same planets before and after
  stripping, chosen by `envelope_loss_fraction`.

### 4.3 Hycean worlds
Madhusudhan, Piette and Constantinou (2021) [S]: an ocean under an H2-rich
atmosphere; radii up to 2.6 Earth radii at 10 Earth masses and 2.3 at 5;
ocean-surface limit about 400 K; inner edge at T_eq up to about 500 K for late
M dwarfs, no outer limit. K2-18 b (8.6 Me, 2.6 Re, T_eq about 255 K): the 2023
CH4 and CO2 detection is established, the 2025 DMS claim is not confirmed, and
Hycean against gas-rich sub-Neptune is open [S].

Generator window for a Hycean subtype of class R [D from S]: M 1 to 10 Me, R 1.4
to 2.6 Re (cap 2.3 below 5 Me), density 1.5 to 3 g/cm3, T_eq 150 to 500 K, H2
pressure at the ocean surface 1 to 1e4 bar (Atm says 10 to 100), ocean-surface
T 280 to 400 K. An ocean of 100 km or more has ice VI or VII at the bottom, so
Master's ice-sealed scoring (`L_chem` about 0.05) applies to deep Hyceans; only
shallow ones (under about 100 km,
[activity-magnetism-radiation-hydrosphere.md](activity-magnetism-radiation-hydrosphere.md)
5.2) score high, which resolves the Master table's tension (0.74 against 0.12).
A Hycean keeps primordial H2 at v_esc 15 to 20 km/s below T_exo of about 300 to
500 K and loses little (`f_lost` about 0.01 percent at S = 1 around a Sun-like
star, 0.15 percent around an M dwarf, for 8.6 Me at 2.6 Re [C]).

## 5. Classes

### 5.1 Proposed definitions

| Code | Name | Zones | Radius (Earth) | Mass (Earth) | Atmosphere rule | Notes |
|---|---|---|---|---|---|---|
| R (new) | Sub-Neptune (+ Hycean flag) | e, c, h | 1.6 to 3.2 (cap 4) | 3 to 10 | H2/He envelope of 0.5 to 10 percent of mass if `envelope_loss_fraction` is under 0.5 of the birth envelope, else S | Rogue yes. Moon no. Density 1.2 to 4 g/cm3 |
| S (landed) | Rocky super-Earth | h, e, c, r | 1.2 to 1.8 | 2 to 10 | None only when `r > 30`; trace for 2.5 to 30; otherwise a secondary atmosphere by the retention function | Airless only above T_eq 650 to 1,200 K [C] |
| N (re-scope) | Desiccated hot world (Venus type) | h (not e) | 0.78 to 1.8 (5,000 to 11,500 km) | 0.3 to 10 | `S_runaway < S < S_cs`; CO2 10 to 100+ bar, N2 1 to 4 bar, SO2/H2SO4; T 450 to 900 K | Lifeless, no life chemical. Needs GEN.85 species |
| Z (new) | Lifeless temperate world | e | 0.8 to 1.6 | 0.3 to 6 | Never-metabolised air: dense CO2/N2, stalled CO (M dwarf), or abiotic O2 10 to 100+ bar (young M dwarf, water lost) | Needs the pre-MS rule and GEN.85 |
| U (new) | Icy world (ice dwarf, large icy moon) | c, r | 500 to 3,000 km | -- | `r < 1` and T below the N2 limit gives 1 to 100 micro-bar N2/CH4; else none | Moon yes. Rogue yes |
| W (new) | Small subsurface-ocean body (Enceladus) | c (moon) | 50 to 500 km | -- | Plume exosphere only | Needs tidal heating; rogue no |
| X (new, extend) | Subsurface-ocean world (Europa to super-Earth) | c, e, r | 500 to 11,500 km | up to 10 | Trace O2/H2O if `r > 2.5`; thicker CO2/N2/CH4 if `r < 1` | Radius ceiling raised from 10,000 km. Rogue yes above about 0.1 Me |
| Y (new) | Titan-like | c | 1,500 to 4,000 km | 0.01 to 0.1 | `r < 0.7`, T_surf 80 to 120 K, N2 0.5 to 2 bar + CH4 | Needs S under about 0.02. Moon yes. Rogue no |
| T (retire) | Gas dwarf | c | 15,000 to 55,000 km | 6 to 40 | -- | Merge into R (below 10 Me) and I (above); frees the letter |
| I (adjust) | Ice giant | e, c, r | 3.2 to 8 | 10 to 60 | H/He + ices | Reduces the overlap with R |

Sweep findings for GEN.29 [C, from the section 2 medians]: A and B hold 13 to
17 kPa on 3.4 to 3.7 km/s bodies at S about 4, about 50 times beyond the
shoreline (such bodies should be airless or hold micro-bar volcanic air); N
sits in the ecosphere at S about 1.05, below the runaway threshold (1.107 for
1 Earth mass [S]); K, L and P are at or below the Armstrong limit (full list in
[habitability-index.md](habitability-index.md) 5); P (r about 0.1) is
consistent with the shoreline, so GEN.27 is independent of retention. A rocky
rogue of 10 to 16 Me should be R, not a 17,600 km Class S: an Earth-like 16 Me
planet is 2.1 Re = 13,400 km, and 4.2 g/cm3 at 2.76 Re is a Neptunian density.

### 5.2 Letters, options and PR order
`planet_class` is `String(16)` on planets and moons and `String(4)` on rogues,
so a longer code is storable, but the simplest plan creates no new letter for
GEN.91.

- **Option A (recommended):** re-scope N (hot, lifeless, up to 1.8 Re), extend
  X, Hycean as an R flag, S gets air by rule, T retired. The "S and V in the
  other zones" classes then exist as N (hot) and X or R (cold).
- Option B: two-character codes (for example `NS`, `XS`); every
  `tuning.PLANET_CLASSES` lookup and the `String(4)` rogue column must handle
  them.
- Option C: retire T and migrate Q; Q cannot be reused safely, so this frees
  one letter.

PR order under GEN.33 (one class per PR):

0. `planetgen/physics/atmosphere_retention.py`: the 3.6 functions with unit
   tests on the 13-body table; no class changes (GEN.91 part 1).
1. R (sub-Neptune), with the envelope-survival rule and the retirement of T.
2. S fix: the atmospheric rule for non-hot S. It changes existing output, so it
   needs GEN.92's "existing galaxies unchanged until regenerated" guard.
3. U (icy world): data only.
4. X (and W): needs the subsurface-ocean data (`rogue_surface.py` has the shell
   and ocean depth for rogues).
5. Y: needs the vapour-pressure gate (3.3).
6. N re-scope: needs the runaway thresholds and removal of N's life chemical.
7. Z: last, because it needs GEN.85 partial pressures and the pre-main-
   sequence rule.
8. Hycean flag on R; 9. GEN.29 sweep.

## 6. Life-stage inputs missing from the code
- **Time since habitable, not star age.** A planet is not habitable before the
  star leaves its pre-main-sequence runaway (M dwarfs: 0.1 to 1 Gyr [R]) and
  before the surface cools (0.05 to 0.1 Gyr): `t_hab = star_age - t_start`.
  X-ray saturation (`t_sat`, 1 to 4 Gyr for M dwarfs) is a different clock;
  Rad merges them (section 8).
- **Tidal locking of planets** [C]: Q = 100, k2 = 0.3 for rock. An Earth-mass
  habitable-zone planet is locked below about 0.45 Msun (middle of the zone,
  10 h initial spin, 5 Gyr) to 0.65 Msun (S = 1, 24 h, 10 Gyr); use the formula,
  not a cut
  ([activity-magnetism-radiation-hydrosphere.md](activity-magnetism-radiation-hydrosphere.md)
  3.3). A locked planet with air under about 0.1 bar risks night-side collapse
  (3.3).
- **Flares and XUV** decide whether air and oceans survive (3.4) and enter the
  score through ozone and `L_rad`.
- **Abiotic O2**: an M-dwarf habitable-zone planet may have lost its water and
  carry O2 that never saw life (Luger and Barnes [S]): Class Z's M-dwarf
  subtype.

## 7. Life stage against time habitable (GEN.92)

### 7.1 Earth against the project's timelines
Earth formed 4.54 Ga. Dates from search results [S: Wikipedia "Timeline of
life"]: earliest evidence of life about 3.8 Ga, oxygenic photosynthesis 3.0 to
3.5 Ga (disputed), Great Oxidation 2.4 to 2.5 Ga, earliest confirmed eukaryote
about 1.85 Ga, first animals about 0.6 Ga. Elapsed time since formation: life
by 0.7 to 1.0 Gyr; oxygenic photosynthesis 1.0 to 1.5; eukaryotes about 2.7;
animals 3.9 to 4.0; technology about 4.5.
`EVOLUTIONARY_TIMELINES["normal"]` (abiogenesis 0.5, photosynthesis 1.5,
complex cells 2.5, multicellularity 4.0, technology 4.5 Gyr) fits closely, for
one sample. The `slow` row (3 to 35 Gyr) and `fast` row (0.005 to 0.08 Gyr) are
invention: Earth bounds only the first stage from above.

`tuning.STAR_EVOLUTION` lifespans follow `t_MS = 10 Gyr M^-2.5` (F 4.3 to 9.1,
G 9.1 to 17.5, K 17.5 to 73.5, M 73.5 to 5,500 Gyr) [C]. O and B stars (0.001
to 1.5 Gyr) give life nothing, yet the repo gives them "Retinal" with fast
evolution.

### 7.2 Rule table
`t_hab` is star age minus the habitable-onset delay (0.05 Gyr; for stars under
0.5 Msun add 0.1 to 1 Gyr of pre-MS runaway if the planet would have been above
`S_runaway` then). Pace multiplier from the existing scales (normal 1, slow 2.5
to 5, fast 0.01 to 0.1 if Boss keeps it). Gates use the scores of
[habitability-index.md](habitability-index.md).

| Stage | Earth-based minimum `t_hab` (normal) | Required gates | Capped by |
|---|---|---|---|
| 0 Sterile | -- | `PHI_bio < 0.15` (no liquid solvent, ice-sealed, dose beyond 10 Gy/yr, T above 122 C) or `t_hab` below the abiogenesis time | Classes N, U with no ocean, Z |
| 1 Abiogenesis (microbial) | 0.5 Gyr | `PHI_bio >= 0.3`; `L_solv >= 0.3`; 258 to 395 K somewhere | Star lifespan >= 0.5 Gyr (excludes O, B) |
| 2 Photosynthesis | 1.5 Gyr | `PHI_bio >= 0.3`; light (`L_ener`) | Star lifespan >= 1.5 Gyr |
| 3 Complex cells | 2.5 Gyr | Stage 2; free O2 about 1 percent of present (0.2 kPa) from oxygenic photosynthesis [R]; anoxygenic-only M-dwarf biospheres stall at 2 unless an oxygenic variant is drawn | Star lifespan >= 2.5 Gyr (M <= about 1.6 Msun) |
| 4 Multicellular | 4.0 Gyr | `PHI_cpx >= 0.3`; pO2 >= 8 kPa; pCO2 < 2 kPa; pCO < 0.01 kPa; 273 to 330 K; dose under 10 Gy/yr | Star lifespan >= 4.0 Gyr (M <= about 1.35 Msun) |
| 5 Technological | 4.5 Gyr | Stage 4; land fraction >= about 5 percent; pO2 >= about 15 kPa for combustion [R]; stable obliquity and climate | Star lifespan >= 4.6 Gyr; optionally an M-dwarf cap at stage 4 while age < `t_sat` |

Highest stage = the latest stage whose `t_hab` is reached on the planet's pace
and whose gate is passed. Today only the time part exists, and the gate is a
per-class cap. Replacing the cap by the gates settles GEN.90's known conflicts.


## 8. Corrections to the source documents
"Atm" = Atmospheric Toxicity.md, "Math" = Mathematical and Algorithmic
Implementation of the Planetary Habitability Index.md, "Master" = Planetary
Habitability and Speculative Xenobiology.md, "Rad" = Naturally Occuring
Ionizing Radiation.md. Threshold conflicts C1 to C16 (mask band, CO2, CO, O2,
wet-bulb, pressure) are resolved in [habitability-index.md](habitability-index.md)
3; the ones below concern atmospheres, classes and timing.

| # | Source and section | Statement | Correction |
|---|---|---|---|
| A1 | Atm "Climatic and Density Extremes", Sub-Baric Ice Worlds; Master archetypes | Low-mass bodies (0.3 to 0.5 Me) near the middle or outer zone lose air by kinetic escape to 10 to 20 kPa | A 0.3 to 0.5 Me rocky planet has v_esc 7.2 to 8.7 km/s and `S_cs` 1.0 to 2.2 against zone insolation 0.3 to 0.7, so r is 0.2 to 0.6 and it retains air. Low pressure needs impact erosion, low supply or cold-trapping [C] |
| A2 | Atm Sub-Baric Ice Worlds; Master | Pure-O2 mask at 10 to 20 kPa "matches alveolar oxygenation" | Alveolar pO2 of pure O2 at 10 to 20 kPa total is -3 to +7 kPa; a mask needs 19 to 21 kPa or more, below that a pressure suit ([habitability-index.md](habitability-index.md) 4) |
| A3 | Atm "Post-Runaway Desiccated Super-Earths"; Master (tech 0.81) | A 1-bar dry "dune world" with masks is the desiccated outcome | Venus type (9 MPa, 737 K) is the default; the dune world is a rare variant (4.2) |
| A4 | Atm Archetype C; Master Hycean row | Hycean 10 to 100 bar H2, `PHI_bio` 0.74 | 1 to 1e4 bar; deep oceans (100 km or more) are ice-sealed with `L_chem` about 0.05 (4.3) |
| A5 | Atm "Hydrodynamic and Kinetic Escape", Archetype B; Master | Abiotic hyperbaric O2 worlds of 10 to 100+ bar at 240 to 300 K | Supported only for habitable-zone planets of M dwarfs older than about 1 Gyr that lost several oceans [S]; make it Class Z's M-dwarf subtype |
| A6 | Atm Archetype A and "HZCL mismatch"; Master | Dense CO2 of 2 to 10 bar; complex-life limit "0.005 to 0.05 bar" | Fine for microbes; the metazoan gate is 0.5 to 5 kPa, so dense CO2 is lethal to metazoans ([habitability-index.md](habitability-index.md) 3, C3, C11) |
| A7 | Math 2.1 escape and Jeans | `Mdot = eps pi R_p R_XUV^2 F / (G M K_tide)`; `lambda_esc,i` with no threshold | Valid; no values given: use `eps` 0.1, `R_XUV` = R_p, `K_tide` 1 (the 3.4 form, a lower bound). Jeans: kept if 30 or more, lost if under 15 (3.2) |
| A8 | Math 4 pseudocode | `P_0 = max((initial_volatiles - atm_loss) * gravity, 0)` | Dimensionally a force: surface pressure is `M_atm g / (4 pi R^2)`. No species are lost, so use the 3.6 function per gas |
| A9 | Math 4 pseudocode, 2.1 | `seafloor_pressure > 1.2e9` returns 0 (ice VI/VII seal) | The seal depends on temperature: 0.2 to 0.35 GPa near 251 to 256 K, 0.63 GPa at 273 K, 2.2 GPa at 355 K [S]. Use the liquidus at the bottom temperature ([activity-magnetism-radiation-hydrosphere.md](activity-magnetism-radiation-hydrosphere.md) 5.2) |
| A10 | Rad "Host Star Spectral Regimes" | M dwarfs undergo "prolonged pre-main-sequence contractions of 1.0 to 2.5 Gyr" with saturated XUV | Conflates PMS contraction (0.1 to 1 Gyr [R]) with X-ray saturation (`t_sat` 1 to 4 Gyr for 0.1 to 0.45 Msun) |
| A11 | PHI-4 / TODO GEN.90 | N "carries life with no stage cap" | N has a life chemical and a hidden uncapped timeline, but is not in `HABITABLE_PLANET_CLASSES` (section 2) |

## 9. Evidence notes
[R] claims to verify when paper access is allowed:

1. Kopparapu `c`, `d` coefficients, the 0.1 and 5 Earth-mass runaway sets and
   the recent-Venus and early-Mars rows (they reproduce the quoted solar
   distances 0.75, 0.95, 0.993, 1.676, 1.774 AU).
2. Zahnle and Catling's shoreline constant (K derived here from 13 bodies);
   Jeans thresholds (30 kept, 15 lost); the exobase model (five bodies).
3. The XUV history in 3.4 (see the activity document's evidence notes).
4. Radius-valley slope, Zeng et al. `R = M^0.27`, the 1.12 factor.
5. Hycean pressure range and the deep-ocean ice-VII argument.
6. Tidal-locking constants (Q 100, k2 0.3), M-dwarf pre-MS runaway durations,
   Joshi et al. 1997's 0.1 bar threshold.
7. LHS 3844 b "bare rock" and the planet parameters used for TRAPPIST-1 and
   55 Cnc e.

A Cockell "nine habitability factors" list and a Barnes habitability index
were searched for and not found; neither is used.

## 10. Sources
Found by search (the full papers could not be opened):

- Shoreline: https://arxiv.org/abs/1702.03386, https://arxiv.org/pdf/2508.12865
  (also https://arxiv.org/pdf/2504.19872, https://arxiv.org/pdf/2507.02136)
- Habitable-zone limits: https://alphaxiv.org/abs/1301.6674,
  https://astrobiology.com/2013/06/26/how-close-is-earth-to-a-runaway-greenhouse/,
  https://arxiv.org/pdf/1910.07573, https://arxiv.org/pdf/1608.06772,
  https://arxiv.org/pdf/1707.07986
- Hycean and K2-18 b: https://arxiv.org/pdf/2108.10888,
  https://arxiv.org/html/2505.13407v1, https://arxiv.org/pdf/2508.05961,
  https://www.aanda.org/10.1051/0004-6361/202555580
- Sub-Neptunes: https://arxiv.org/abs/1603.08614 (code and parameters read from
  https://raw.githubusercontent.com/chenjj2/forecaster/master/),
  https://arxiv.org/pdf/1710.05398, https://sseh.uchicago.edu/doc/Owen_and_Wu_2017.pdf,
  arXiv:1407.4457
- Water loss and observed atmospheres: https://arxiv.org/abs/1411.7412,
  https://arxiv.org/abs/2405.04744, https://www.arxiv.org/pdf/2412.11627,
  https://mpg.de/23861226/ducrot_trappist-1b_natureastronomy_2024.pdf,
  https://astrobiology.com/2025/09/jwst-tst-dreams-nirspec-prism-transmission-spectroscopy-of-the-habitable-zone-planet-trappist-1-e.html
- Earth timeline: https://en.wikipedia.org/wiki/Timeline_of_life
- Pressure breathing sources are listed in [habitability-index.md](habitability-index.md) 7.2 and 7.3.

