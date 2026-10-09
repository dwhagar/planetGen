# Stellar activity, magnetic fields, surface radiation dose and hydrosphere

The physics inputs the habitability index needs beyond atmosphere species:
how active each star is, whether a planet has a magnetic field, how much
ionizing radiation and UV reaches its surface, and what water and ocean it
holds. This is the research behind `docs/TODO.md` items GEN.86, GEN.87 and
GEN.88, the spin dependency GEN.104, and the reuse of
[rogue-planet-surface.md](rogue-planet-surface.md). It carries the still-valid
formulas of Boss's uploaded documents "Naturally Occuring Ionizing
Radiation.md" (Rad), "Chemical Habitability.md" (Chem), "Planetary
Habitability Index.md" (PHI-4) and "Mathematical and Algorithmic
Implementation of the Planetary Habitability Index.md" (Math), and names each
correction to them (section 6). Atmosphere retention and class definitions are
in [atmospheres-retention-and-classes.md](atmospheres-retention-and-classes.md);
the score structure and tiers are in [habitability-index.md](habitability-index.md).

Informs: GEN.84, GEN.86, GEN.87, GEN.88, GEN.89, GEN.104, rogue-planet-surface.md
Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [S] seen in a search result (URL under Sources), [C] computed
(inputs shown), [R] recalled from memory or a design default, unconfirmed.
The research environment could only read search-result text, not the papers;
the [R] items are listed under Evidence notes.

## Decisions already taken

- Boss (2026-10-03 05:38Z): "Create a habitability index based on the
  pressure, temperature, composition, etc..." (GEN.83).
- Boss (2026-10-07 17:11Z, confirmed 2026-10-09 07:48Z): PHI-4's four domains
  and colour tiers for display, with the Xenobiology doc's three scores
  (`PHI_bio`, `PHI_cpx`, `Phi_tech`) behind them
  ([habitability-index.md](habitability-index.md)). The radiation, magnetism
  and ocean numbers below feed those scores.

## 1. Summary

1. **Stellar activity needs four numbers per star** (2.1): saturated
   `L_X/L_bol` (7.4e-4), saturation age `t_sat`, a decay exponent and a flare
   rate at 1e33 erg. The Rad doc's "L_XUV/L_bol about 1e-3" is the X-ray band
   only; total XUV of a saturated star is 3.3e-3 (G) to 7.3e-3 (0.1 Msun) [C].
   Habitable-zone planets of M dwarfs get 1,000 to 2,000 times Earth's present
   XUV for billions of years; an exposure index gives the tiers (2.3).
2. **Magnetic fields** follow core size, density and convective power, nearly
   independent of rotation while the dynamo stays dipolar. A tidally locked
   M-dwarf planet is probably multipolar (dipole about 0.15, standoff 0.53x).
   **Code finding:** only moons are tidally locked in
   `generate_orbital_motion_properties`; an M-dwarf-zone planet gets a random
   10 to 1,400 h spin but is locked within 15 kyr (0.1 Msun) to 3.4 Myr
   (0.3 Msun) [C]. GEN.104 must fix this before GEN.86.
3. **Surface dose** is dominated by column mass. Mars surface 0.64 mSv/day dose
   equivalent; Earth's cosmic dose 0.39 mSv/yr, 2.4 being the total natural
   background [S]. A fitted calculator (4.2) matches within 10 percent.
4. **Hydrosphere:** the rogue-planet ice-shell function reproduces Europa and
   Enceladus [C], so reuse is sound with four fixes (5.3). The main one is a
   **high-pressure-ice cap**: a liquid ocean cannot exceed about 70 to 110 km at
   Earth gravity, so `ocean_depth_km` up to 855 km is not liquid to the bottom.
5. Add a UV/ozone flag and an optional galactic-hazard modifier to GEN.87 (4.5).

## 2. Stellar activity (GEN.86)

### 2.1 X-ray and XUV evolution

- **Saturation level.** Wright et al. 2011 (824 stars): `log(L_X/L_bol)` =
  -3.13 +/- 0.08, saturation at Rossby number 0.13, unsaturated slope -2.70
  [S for slope and Rossby number; -3.13 is [R]]. Jackson, Davis and Wheatley
  2012 (717 stars): the saturated ratio falls from about 10^-3.1 (late K) to
  10^-4.3 (early F) [S], so "10^-3 for all stars" holds only for K and M.
- **M dwarfs saturate longer.** West et al. 2008: H-alpha active lifetimes rise
  from about 0.8 Gyr (M0) to about 8 Gyr (M7) [S]; an upper bound for X-ray
  saturation.
- **Spread at fixed age is real.** Tu et al. 2015: a solar-mass star stays
  saturated about 10 Myr if it starts slow, 300 Myr if fast [S]; draw a
  per-star offset of 0.3 to 0.5 dex.
- **Decay.** Skumanich spin-down gives `L_X ~ t^-1.35 to -1.5` [C], too shallow
  for the Sun today. Exponent 2.0 for M >= 0.5 Msun gives 3.5e-7 for the Sun now
  and 2.1e-5 at 0.6 Gyr (Hyades about 1e-5 [R]); fully convective stars decline
  slowly (Wright 2018 [R]), so use 1.0 below 0.3 Msun.
- **XUV from X-ray.** `log L_EUV(100 to 920 A) = 4.80 + 0.860 log L_X` (erg/s)
  [R: Sanz-Forcada 2011 as used by Chadney 2015]. It gives `L_XUV(Sun now)` =
  1.5e28 erg/s, 5.3 erg/cm2/s at 1 AU against Ribas's 4.6 [R], and
  `L_XUV ~ t^-1.29`, matching the -1.2 to -1.3 band exponents of Ribas et al.
  2005 [S].

Defaults [C: design defaults calibrated as above]:

| Mass (Msun) | Type | `L_X/L_bol` saturated | `t_sat` (Gyr) | Decay exponent | `L_XUV/L_bol` saturated |
|---|---|---|---|---|---|
| 1.0 to 1.4 | F, G | 7.4e-4 (G); falls to 5e-5 by 1.4 | 0.10 | 2.0 | 3.3e-3 |
| 0.8 | K | 7.4e-4 | 0.15 | 2.0 | 3.6e-3 |
| 0.6 | K/M | 7.4e-4 | 0.4 | 2.0 | 4.1e-3 |
| 0.45 | M | 7.4e-4 | 1.0 | 1.75 | 5.0e-3 |
| 0.3 | M | 7.4e-4 | 2.0 | 1.0 | 5.6e-3 |
| 0.2 | M | 7.4e-4 | 3.0 | 1.0 | 6.0e-3 |
| 0.1 | late M | 7.4e-4 | 4.0 | 1.0 | 7.3e-3 |

Validation: the Sun at 4.57 Gyr gives `L_X` 1.4e27 (target 1e27 to 2e27);
Proxima (0.12 Msun) 3.8e27 against about 4e26 to 1.7e27 observed [R], so late M
dwarfs are about 2x too active. Add a log-normal per-star offset of 0.4 dex
after `t_sat`. This model replaces the stand-in XUV history in
[atmospheres-retention-and-classes.md](atmospheres-retention-and-classes.md)
3.4.

### 2.2 XUV flux at the planet

`F_XUV = L_XUV / (4 pi d^2)`. With the project's mass-luminosity law
(`tuning.MS_MASS_LUMINOSITY_PIECES`) and a planet at S = 1 (d = sqrt(L) AU),
in units of Earth's present 4.6 erg/cm2/s [C]:

| Star (Msun) | d (AU) | 0.1 Gyr | 0.5 Gyr | 1 Gyr | 2 Gyr | 4.5 Gyr | 8 Gyr |
|---|---|---|---|---|---|---|---|
| G 1.0 | 1.00 | 974 | 56 | 17 | 4.9 | 1.2 | 0.4 |
| K 0.8 | 0.64 | 1,074 | 128 | 38 | 11 | 2.7 | 1.0 |
| K/M 0.6 | 0.36 | 1,224 | 825 | 243 | 72 | 17 | 6.4 |
| M 0.4 | 0.167 | 1,465 | 1,465 | 1,465 | 765 | 263 | 124 |
| M 0.2 | 0.075 | 1,776 | 1,776 | 1,776 | 1,776 | 1,245 | 752 |
| M 0.1 | 0.034 | 2,165 | 2,165 | 2,165 | 2,165 | 1,953 | 1,182 |

The saturated flux at S = 1 is nearly constant; M-dwarf planets differ in how
long they stay there. Scale by S otherwise.

### 2.3 Effects, exposure index and ozone flag

- **Erosion.** Energy-limited loss `Mdot = eps pi F R^3 / (G M K)`. With
  eps = 0.1 and 6.26e7 J/kg to lift water, an Earth-size planet at S = 1 loses
  4.5 ocean-equivalents by 5 Gyr, yet Earth kept its oceans because the cold
  trap limits real loss [C]. Use it as a relative stress index, and as a real
  estimate only past the moist-greenhouse limit. Scale by
  `(R/R_E)^3 / (M/M_E)`.
- **Ozone.** Segura et al. 2010: with protons, a large flare removes 94 percent
  of ozone over two years on an unmagnetized planet, recovering over about
  50 years; Howard et al. 2018: 90 percent loss in five years for repeated
  Proxima-type flaring [S].
- **Ozone-loss flag** [R, design default]: flare irradiation index
  `FI = N33 / d_AU^2` (Earth 0.003; Proxima b 2,100). Flag when `FI * f_B >= 300`,
  `f_B` = 1 (no dipole), 0.3 (magnetopause under 3 R_p), 0.1 (above), so that
  Proxima b without a field is flagged and Earth is not.

Exposure index: cumulative XUV to age t relative to a G star at 1 AU over
5 Gyr (= 4.5 ocean-equivalents) [C]:

| Star | to 0.5 Gyr | 1 | 2 | 5 | 10 |
|---|---|---|---|---|---|
| G 1.0 | 0.9 | 0.9 | 1.0 | 1.0 | 1.0 |
| K 0.8 | 1.3 | 1.5 | 1.6 | 1.7 | 1.7 |
| K/M 0.6 | 2.8 | 3.8 | 4.4 | 4.9 | 5.1 |
| M 0.45 | 3.3 | 6.6 | 10.5 | 13.8 | 15.4 |
| M 0.3 | 3.7 | 7.5 | 15.0 | 29.6 | 41.8 |
| M 0.1 | 5.0 | 10.2 | 20.5 | 50.4 | 81.1 |

Tiers: 2 or less Earth-like, 2 to 8 elevated, 8 to 30 high, over 30 extreme.
The shoreline research multiplies this index by S in place of bolometric S.

### 2.4 Flares

- **Sun-like stars.** Maehara et al. 2012: one flare above 1e34 erg per 800
  years on solar types, mostly faster rotators than the Sun [S]. Notsu 2019:
  slow rotators once per 2,000 to 3,000 years [S]. Okamoto 2021: 1e34 every
  6,000 years for the Sun [S].
- **Proxima.** 1e33 erg about 8 per year (MOST, extrapolated), at least 5.2
  (Evryscope), 3 (TESS) [S]; use 3 to 8.
- **Slope** of the frequency distribution 1.4 to 2.4, 1.99 +/- 0.07 for
  low-activity M dwarfs, no change with age (Davenport 2019) [S]. Draw alpha
  uniformly in 1.8 to 2.2. **Flaring fractions** (TESS, Gunther 2020): about
  30 percent of mid to late M, 5 percent of early M, under 1 percent of F, G, K
  [S].

Generator table, `log10 N(>1e33 erg)` per year [R/C]. The G cell (2 to 6 Gyr)
of -2.5 gives 1e34 once per 3,000 years at alpha = 2 (Notsu, Okamoto); the
late-M cell (2 to 6 Gyr) of 0.8 contains Proxima (0.7); the rest interpolate,
+/- 0.5 dex.

| Class (mass, Msun) | <0.1 Gyr | 0.1 to 0.6 | 0.6 to 2 | 2 to 6 | >6 |
|---|---|---|---|---|---|
| G (0.85 to 1.1) | 1.0 | 0.0 | -1.3 | -2.5 | -2.8 |
| K (0.6 to 0.85) | 1.3 | 0.5 | -0.8 | -2.0 | -2.3 |
| Early M (0.35 to 0.6) | 2.0 | 1.6 | 0.8 | -0.3 | -0.7 |
| Mid M (0.15 to 0.35) | 2.3 | 2.2 | 1.8 | 0.9 | 0.5 |
| Late M (<0.15) | 2.5 | 2.5 | 2.2 | 0.8 | 0.5 |

A, F and B stars get `log N33` = -4. Flares above 1e34 erg subtract
`(alpha - 1)` dex. Flares scale the particle dose (4.2).

## 3. Planetary magnetic fields (GEN.86)

### 3.1 The scaling law

Olson and Christensen 2006 fitted 145 numerical dynamos; the time-averaged
dipole moment scales with buoyancy flux as `F^(1/3)` in the dipolar regime and
predicts solar-system moments to about 15 percent [S via a citing paper]:

```
M_dip = 4 pi r_o^3 gamma sqrt(rho / (2 mu0)) (F D)^(1/3),   gamma about 0.2   [R: recalled form]
```

(`r_o` outer-core radius, `rho` core density, `D` shell thickness, `F` buoyancy
flux in m2/s3.) Check [C]: `r_o` 3480 km, `rho` 11,000 kg/m3, `D` 2260 km,
`F D` = 1.4e-6 gives 7.8e22 A m2 (Earth).

Rotation: the Christensen and Aubert 2006 law `Lo = 0.92 sqrt(f_ohm) Ra_Q*^0.34`
[S] with `Ra_Q* ~ Omega^-3` gives `B ~ Omega^-0.02`, so field strength does not
depend on rotation in the dipolar regime (a derivation here). Rotation sets the
local Rossby number `Ro_l`: dipolar below about 0.12, multipolar above, with
`Ro_l ~ Omega^-1.23` [R]. A hard limit would fall at P about 30 hours for
Earth, too sharp (Mercury, 59 d, is dipole-ish), so use a soft switch (3.5).

### 3.2 Measured moments

Dipole moments in A m2 (`M = 4 pi R^3 B_eq / mu0`, [C from recalled `B_eq`]):
Earth 7.8e22 | Mercury 2.8e19 | Ganymede 1.3e20 | Jupiter 1.56e27 | Saturn
4.6e25 | Uranus 3.8e24 | Neptune 2.1e24. Venus and Mars have no active dynamo;
Io, Europa, Callisto and Titan have none of their own (Europa and Callisto have
induced fields from their oceans, a flag for GEN.88).

### 3.3 Tidal locking

Gladman et al. 1996 (the form used by `_tidal_locking_timescale_seconds`):
`t = (2Q / 15 k2) omega0 a^6 m / (G M*^2 R^3)`. The project uses Q/k2 = 3,333
for moons (Q 100, k2 0.03). For a rocky planet use k2 = 0.3 (10x faster).
Earth-mass, Earth-size planet at S = 1, initial P = 24 h [C]:

| Star (Msun) | a (AU) | `t_lock`, k2 = 0.3 | `t_lock`, project Q/k2 |
|---|---|---|---|
| 1.0 | 1.00 | 1.0e11 yr | 1.0e12 |
| 0.8 | 0.64 | 1.1e10 | 1.1e11 |
| 0.6 | 0.36 | 6.1e8 | 6.1e9 |
| 0.45 | 0.20 | 3.4e7 | 3.4e8 |
| 0.3 | 0.12 | 3.4e6 | 3.4e7 |
| 0.2 | 0.075 | 4.6e5 | 4.6e6 |
| 0.1 | 0.034 | 1.5e4 | 1.5e5 |

Habitable-zone planets are locked for stars below about 0.6 to 0.7 Msun at
these inputs (nearer 0.45 Msun with a 10 h initial spin and the middle of the
zone); use the formula, not a mass cut. The project's formula matches Gladman
(the Observational Kinetics doc's 4/9 prefactor is 4/3 higher; immaterial).

### 3.4 Magnetosphere standoff

Pressure balance with a dipole: `R_mp = [mu0 f0^2 M^2 / (8 pi^2 p_sw)]^(1/6)`,
`f0` = 1.16 (other papers 1.3), `p_sw` the total wind pressure [S]. For Earth
(M 7.768e22, p_sw 2.24 nPa) this gives 9.76 R_E [C], matching the empirical
`R/R_E = 9.75 (M/M_E)^(1/3) (P/P0)^(-1/6)`. A tenfold wind pressure shrinks
`R_mp` only 1.47x [S].

`p_sw = 2.24 nPa * (Mdot per area / solar) / d_AU^2`, with `Mdot/area ~ F_X^1.3`
(Wood-type [R]), capped at 10x solar for G and K and 1x for M under 0.4 Msun
(Proxima's is about 0.2 solar per area [R]). Earth-dipole planet at S = 1 [C]:

| Star (Msun) | `p_sw` / Earth's | `R_mp` (R_p) |
|---|---|---|
| 1.0, young / 4.5 Gyr | 10x / 1x | 6.6 / 9.7 |
| 0.6 | 54x | 5.0 |
| 0.4 | 36x | 5.4 |
| 0.2 | 176x | 4.1 |
| 0.1 | 870x | 3.2 |

With the multipolar factor 0.15 the radius drops by `0.15^(1/3)` = 0.53, to 1.7
to 3.5 R_p for locked planets, matching Zuluaga 2013's 1.5 to 4.0 (locked) and
3 to 8 (1-day rotators) [S]. CMEs raise `p_sw` by 10 to 1,000 transiently [R].

### 3.5 Generator rules (design defaults [R])

For a rocky planet with mass M (Earth masses), age, rotation period P, core
mass fraction c (Earth 0.325) and lid regime:

1. **Eligibility.** M >= 0.05 and c >= 0.15, else no intrinsic field (Mercury,
   0.055 Me, has one because its core is 70 percent of its mass).
2. **Dynamo lifetime.** `tau_dyn = 4.0 Gyr * M^0.5 * U(0.5, 1.5)`; times 0.4 for a
   stagnant lid (Venus, Mars) and 0.5 if c < 0.2. Active if age < `tau_dyn`.
   Zuluaga and Cuartas-Restrepo: low-mass super-Earths have strong fields
   lasting 2 to 4 Gyr when P > 1.5 d; the literature is unsettled [S].
3. **Strength.** `M_dip = 7.8e22 * M^1.0 * (c/0.325)^1.0 * LN(0.4 dex) *
   (1 - age/tau_dyn)^0.3`. Mercury-like (M under 0.1, c over 0.5): 3e19 *
   LN(0.5 dex).
4. **Rotation switch.** `Ro_l = 0.09 (P / 24 h)^1.23`; `f_dip = 0.15 + 0.85 /
   (1 + (Ro_l / 0.12)^4)`. A locked M-dwarf planet has `f_dip` about 0.15.
5. **Giants.** Always magnetized; log-normal 0.4 dex by bin: Neptune and
   sub-Neptune 2e24 to 4e24; Saturn-mass 5e25; Jupiter-mass 1.5e27, rising
   about linearly with mass above that [R].
6. **Moons.** Ganymede-like (mass >= 0.02 Me, differentiated, resonant tidal
   heating): P = 0.1, `M_dip` = 1.3e20 * LN(0.7 dex). Otherwise none; flag an
   induced field if a subsurface ocean exists.
7. **Output fields:** `magnetic_moment_a_m2`, `dipole_class` (none, weak
   under 0.01 Earth, Earth-like, strong over 5 Earth, multipolar),
   `magnetopause_rp`.

## 4. Surface radiation dose (GEN.87)

### 4.1 Anchors

| Place | Measured | Source |
|---|---|---|
| Earth, sea level, cosmic | 0.39 mSv/yr | UNSCEAR 2008 [S] |
| Earth, all natural | 2.4 mSv/yr (radon 1.26, cosmic 0.39, ground gamma 0.48, ingestion 0.29); about 3.0 in newer UNSCEAR | UNSCEAR [S] |
| ISS | 0.731 mSv/day (GCR only 0.523) | Chang'E-4 LND comparison [S] |
| Moon surface | 1.369 mSv/day (about 500 mSv/yr) | Chang'E-4 LND [S] |
| Mars surface | 0.64 +/- 0.12 mSv/day dose equivalent; 0.21 +/- 0.04 mGy/day absorbed | MSL RAD, Hassler 2014 [S] |
| Europa surface | about 5.4 Sv/day | [R]; searches found only depth profiles [S] |

Mars is 230 mSv/yr (dose equivalent) = 0.077 Gy/yr (absorbed). Compare dose
equivalent with human limits; the `L_rad` term uses absorbed dose in Gy.

### 4.2 The calculator

1. **GCR against column X (g/cm2), no field, 2-pi surface:** log-log
   interpolation through Moon (X 0, 1.37 mSv/day), Mars (16.4, 0.65), a
   high-altitude point (100, 0.30 [R]), then exponential with attenuation
   length 164 g/cm2 [R], fitted so sea level (1033) gives 0.37 mSv/yr against
   0.39 observed. Floor the tail at 1e-3 of the sea-level value for X up to 1e5,
   since muons penetrate kilometres of rock [R].
2. **Magnetic effect on GCR:** global-average factor `1 - 0.4 * min(1,
   log10(Rc0 / 1 GV) / log10(15))`, with cutoff rigidity `Rc0 = 14.9 GV *
   (M_dip / 7.768e22) * (R_E / R_p)^3` [R]. An Earth dipole cuts GCR dose to
   0.6, consistent with Atri 2017 and Griessmeier 2015 (removing the dipole
   raises dose under 2x at 1 bar) [S via Rad]. The "two or three orders of
   magnitude" Rad gives applies to solar particles, not GCR.
3. **Solar and stellar particles (SEP):** free-space annual dose = `20 mSv/yr *
   (N33 / N33_sun)^0.5 / d_AU^2` [R, design]; the square root keeps Proxima b
   from receiving hundreds of Sv/yr. Surface dose = free dose
   `* exp(-p_cut(X) / R0) * polar-cap fraction`. `p_cut` is the rigidity of the
   proton whose range equals X (`R = 0.0022 E^1.77 g/cm2`, E in MeV, within
   2 percent of Bethe-Bloch to 200 MeV [C]); `R0` = 0.2 GV is calibrated on the
   1956 ground-level event [R]. Polar-cap fraction
   `1 - sqrt(1 - sqrt(R_eff / Rc0))`, `R_eff = max(0.5 p_cut, 0.2 GV)` [R].
4. **Crust and radon:** `D_crust = 0.48 mSv/yr * a_rad * H(age) / H(4.5 Gyr)`,
   with `H` the radiogenic heating from `tuning.ROGUE_RADIOGENIC_ISOTOPES`
   (relative 3.5 at 0.5 Gyr, 1.9 at 2, 1.0 at 4.5, 0.51 at 10 [C]). Radon adds
   1.26 mSv/yr (Earth) only with air and porous rock; 0.5 outdoors.
5. **Shielding by ground or water:** the GCR(X) curve with `X = rho * depth`:
   1 m of regolith (1.5 g/cc) is X = 150, 10 m of water X = 1000 [R]. Muons set
   the floor.

### 4.3 Heliosphere compression (GEN.75 link)

The GEN.87 text raises GCR dose in a compressed heliosphere. The effect is real
but small: about 2 to 3x in free-space dose equivalent, and none above X of
about 200 g/cm2, where primaries above 1 GeV dominate [R]. Recommended: `helio_mult = 1 + 1.5 * exp(-X / 100) * min(1, compression)`,
where compression runs 0 to 1 as the heliosphere shrinks from its open-space
radius to under 1 AU. `compressed_heliosphere_radius` in `generation/star.py`
already returns the radius; use R / R_open.

### 4.4 Calculator results (mSv/yr)

| Case | GCR | SEP | Total |
|---|---|---|---|
| Earth sea level (crust and radon not added) | 0.22 | 2e-5 | 0.22 |
| Mars surface, no field | 235 | 0.5 | 236 (0.65 mSv/day, matches MSL) |
| Moon surface | 500 | 20 | 520 |
| 0.05 bar, Earth gravity, Earth dipole | 88 | 0.03 | 88 (Rad: 50 to 100) |
| Proxima-b-like, 1 bar, dipole x0.15 | 0.33 | 1.3 | 1.6 |
| Proxima-b-like, 1 bar, no dipole | 0.37 | 2.6 | 3.0 |
| Proxima-b-like, 0.1 bar, no dipole | 108 | 2,300 | 2,400 (2.4 Sv/yr) |

The 1-bar rows agree with Rad that a dense atmosphere decouples the surface
from flares (under 5 mSv/yr). The 0.1-bar row shows how fast it breaks: a thin
atmosphere around an active M dwarf with no dipole receives Sv per year.

### 4.5 Galactic hazards and UV

- **Supernovae.** Gehrels et al. 2003: a core-collapse supernova within about
  8 pc roughly doubles biologically active UV by ozone loss; the rate within
  8 pc at the Sun is about 1.5 per Gyr [S]. Rule: lethal events per Gyr within
  radius r = `1.5 * (r / 8 pc)^3 * (Sigma_SFR(R) / Sigma_SFR(8 kpc))`; near the
  Galactic centre about 250x higher [R/C]. Gowanlock et al. 2011: the inner
  Galaxy sterilizes more planets, yet a habitable planet is still 10x likelier
  there than in the outer Galaxy [S]. Lineweaver et al. 2004: the Galactic
  habitable zone is an annulus 7 to 9 kpc wide [S].
- **Gamma-ray bursts.** Thomas et al. 2005: a burst within 2 kpc (100 kJ/m2 over
  10 s) destroys 35 percent of ozone globally; a 50 percent cut triples UVB
  [S]. Quoted radii vary (1.4 to 14 kpc); Spinelli et al. 2021: GRBs dominate
  beyond about 2 kpc, supernovae inside [S]. Rule: ozone-loss flag at gamma
  fluence >= 30 to 100 kJ/m2.
- **AGN and quasars.** No source found. A 1e44 erg/s X-ray source on for 1 Myr
  delivers 2.7e5 J/m2 at 3 kpc [C], so a luminous nucleus within a few kpc
  destroys unshielded ozone.
- **UV.** Blackbody 200 to 300 nm flux relative to the Sun: 2.8x (F0), 0.54
  (K0), 0.043 (M0), 0.0011 (M6) [C]; real M dwarfs have chromospheric UV well
  above that, so floor it at 0.05 to 0.2 of the Sun's [R]. Surface DNA-weighted
  UV relative to Earth = `(UV_rel * S) * (O3 / O3_Earth)^-1.6`, clamped near
  1e3 at O3 = 0 (exponent from "a 50 percent ozone cut triples UVB" [S]).

### 4.6 Dose limits and tiers

- Humans: ICRP 20 mSv/yr averaged over five years (occupational), 50 mSv in any
  one year, 1 Sv acute sickness, LD50/30 of 4 to 5 Sv (2.5 to 4.5 Gy depending
  on care) [S for LD50; limits R].
- Microbes: D10 (dose for 90 percent kill) E. coli about 0.3 kGy; spores
  2.1 kGy; Deinococcus radiodurans 10.4 kGy on average [S]; tardigrade LD50 1.1
  to 6.2 kGy [S]. Rad's acute-survival figures are consistent (a different
  metric); use D10 5 to 10 kGy for Deinococcus in any biology table.
- Rad's and PHI-4's bands, written numerically: Blue below 50 mSv/yr, Green 50
  to 100, Yellow "extreme: subsurface", Red above 10 Gy/yr. **There is a gap
  between 100 mSv/yr and 10 Gy/yr** (Mars surface 236 mSv/yr, the Moon 520).
  Proposal for GEN.84 (display bands in
  [habitability-index.md](habitability-index.md) 2): Blue 20 mSv/yr or less
  (ICRP; 50 if Boss prefers the doc's figure), Green up to 100 mSv/yr, Yellow
  (human) 0.1 to 1 Sv/yr, Yellow (microbial) 1 to 10 Gy/yr, Red over 10 Gy/yr.
- The 10 Gy/yr threshold is Eigen's, about the origin of replicating molecules,
  not the survival of extremophiles, which tolerate 1,000x more acute dose.

## 5. Hydrosphere and ocean chemistry (GEN.88)

### 5.1 Inventory, depth and land

Earth's ocean mass is 1.4e21 kg = 2.34e-4 of the planet, mean depth over the
sphere 2.74 km [C]. Mantle water is perhaps 1 to 10 ocean masses and unmeasured
[R].

**Land fraction.** Cowan and Abbot 2014: a tectonically active planet of any
mass keeps continents if its water mass fraction is below about 0.2 percent;
super-Earths store more water in the mantle, whose capacity is the largest
unknown [S]. A toy hypsometry calibrated to Earth's 0.71 ocean fraction (relief
10.5 km at Earth gravity): `f_ocean = min(1, sqrt(2 w M / (rho_w A L)))`,
`L = 10.5 km * (g_E / g)`, `w` the water mass fraction, `A` surface area [C].
It needs a mantle buffer (divide effective surface water by about 4) to reach
the 0.2 percent threshold.

| Water fraction w | Earth gravity | 2.2 g (5 Me, 1.5 Re) | Earth gravity, buffer x4 |
|---|---|---|---|
| 1e-5 | 0.15 | 0.33 | 0.07 |
| 1e-4 | 0.47 | 1.00 | 0.24 |
| 2.3e-4 (Earth) | 0.72 | 1.00 | 0.36 |
| 1e-3 | 1.00 | 1.00 | 0.75 |
| 2e-3 | 1.00 | 1.00 | 1.00 |

Generator: `f_land = 1 - f_ocean` with a buffer log-uniform 1 to 8; a
stagnant-lid planet takes buffer 1. `ROGUE_WATER_MASS_FRACTION_RANGE` (1e-4 to
0.1, log-uniform) makes about 68 percent of water-rich rogues exceed Earth by 4x
and a third exceed 1 percent.

### 5.2 Ocean depth before high-pressure ice

Wagner 2011 triple points, confirmed [S] to 1 to 2 percent: Ih-L-III 251 K at
207.5 MPa; III-L-V 256 K at 346 MPa; V-L-VI 273.3 K at 626 MPa; VI-L-VII 355 K
at 2.2 GPa. With effective water density 1,050 kg/m3 the deepest liquid ocean
(bottom at the liquidus) is [C]:

| g (m/s2) | bottom at 280 K | 300 K | 350 K |
|---|---|---|---|
| 1.3 (Europa) | 553 km | 835 | 1,541 |
| 3.7 | 194 | 294 | 541 |
| 9.8 (Earth) | 73 | 111 | 204 |
| 15 | 48 | 72 | 134 |
| 20 | 36 | 54 | 100 |

An Earth-size ocean of 100 km needs a water mass fraction of 0.85 percent.
Noack et al. 2016: high-pressure ice forms at about 170 km depth on an
Earth-mass planet, but rock heat can melt it from below and layers of nearly
1,000 km are possible if the seafloor is hot [S]. "Deep ocean with rock
contact" and "ice-sealed" are therefore separate outcomes; ice-sealed applies
where bottom pressure exceeds the liquidus.

### 5.3 Reusing the rogue-planet functions

The conductive-lid function `ice_shell_thickness_km(flux, T_top, T_base)`
(`D = (A/F) ln(T_base / T_top)`, A = 567 W/m) reproduces measured shells at
plausible fluxes [C]:

| Body | T_top | Flux used | Result | Reported |
|---|---|---|---|---|
| Europa | 100 K | 0.02 to 0.05 W/m2 | 11 to 28 km | 3 to 30+ km; Juno 29 +/- 10 km [S] |
| Enceladus | 75 K | 0.02 to 0.05 | 15 to 37 km | 20 to 25 km mean, under 5 km at the south pole [S] |
| Ganymede | 110 K | 0.003 to 0.01 | 52 to 172 km | 25 to 150 km [S] |

Tidal heating for generator draws: `E = (21/2) (k2 / Q) G M_p^2 n R^5 e^2 / a^6`
(Peale) [R]. With k2/Q = 0.015 (Io's need) Io gives 9.3e13 W (observed about
9e13) but Europa 7e12 W (estimates 1e11 to 1e12), so use k2/Q = 3e-4 for ice
(Europa 1.4e11 W). Enceladus gives 4.5e8 W against about 1.6e10 W observed, so
allow a resonance boost of 1 to 50x for small icy moons [C]. The radiogenic
flux of a Europa-mass body is 0.0078 W/m2 [C].

Four fixes before sharing the functions with star systems (all verified in
`src/planetgen/physics/rogue_surface.py`):

1. **High-pressure ice cap.** `rogue_surface_conditions` returns `ocean_depth_km`
   up to 855 km (10 percent water). Cap the liquid depth by the 5.2 table and
   report the rest as `hp_ice_km`; for a deep layer the bottom is ice VI/VII and
   the "ocean in contact with rock" assumption (the `has_liquid_water` flag and
   the chemistry classes) is false.
2. **Units.** `ice_shell_thickness_km` is ice thickness, but `water_layer_
   depth_km` is water-equivalent (1,000 kg/m3). The comparison `shell_km >=
   depth_km` and the subtraction overstate the ocean by about 9 percent (ice
   917 kg/m3). Divide the shell by 1.09 when comparing.
3. **Melting at the base.** The base is fixed at 273.15 K. Pressure lowers it by
   0.0074 K/bar (about -10 K for a 100 km shell); salinity and ammonia lower it
   2 to 100 K [R]. `base_temperature_k` is already an argument.
4. **Convection and geometry.** Shells thicker than about 20 to 30 km convect
   and the conductive formula overestimates [R]; on small bodies the lid is a
   large fraction of the radius (Enceladus 20 km of 252 km). Both are within the
   formula's spread.

`water_boiling_k` caps at the critical point (647 K, 22 MPa), which is correct;
above it compressed liquid is stable up to the VI/VII melting curve, which the
function ignores. Move the ice-shell, boiling and water-depth functions into a
shared module (for example `physics/hydrosphere.py`) that both `rogue_surface`
and star-system planets import.

### 5.4 Ocean chemistry classes

Chem defines four classes by pH and water activity. A deterministic assignment
from the GEN.85 state beats a free random draw:

| Class | Conditions | Typical draw |
|---|---|---|
| Ice-sealed | Bottom pressure above the liquidus (depth over the 5.2 table) | pH 3.0 to 5.5, a_w > 0.99 |
| Acid sulfate | Oxidized mantle, SO2 over about 1 percent of volcanic outgassing, low weathering (hot, stagnant lid) or no carbonate sink | pH 1.0 to 4.5, a_w 0.90 to 0.98 |
| Soda | Reduced to intermediate mantle, pCO2 over about 0.1 bar, mafic or ultramafic crust, small continental area | pH 9.0 to 11.5, a_w 0.92 to 0.99 |
| Chloride brine | Water fraction under about 2e-5, or cold and evaporating (endorheic), or cryo-brine margin T under 273 K | pH 3.5 to 6.5, a_w 0.30 to 0.65 |
| Neutral | Everything else (Earth-like: pH 8.1, 35 g/kg, a_w 0.98 [R]) | pH 6.5 to 8.5, salinity 20 to 50 g/kg |

The acid-sulfate range is Chem's (1.0 to 4.0 in the text, 1.0 to 4.5 in the
table); PHI-4's Red band for hyperacidic oceans (pH 0 to 1) is the extreme tail.

Carbonate-silicate buffering (Walker, Hays and Kasting 1981): weathering rises
about as `pCO2^0.3` and about 15 percent per degree [S]. On a stagnant lid it
is supply-limited and can fail, though burial and decarbonation (Foley and
Smye 2018) may restore a feedback with windows of roughly 1 to 5 Gyr [S]. Rule:
a stagnant lid shortens the habitable window to 1 to 5 Gyr and pushes the ocean
to soda or acid sulfate.

### 5.5 Hycean link

Madhusudhan et al. 2021 [S]: massive oceans under H2-rich air, up to 2.6 R_E and
10 Me; inner edge at T_eq up to about 500 K for late M, no outer limit.
`rogue_surface`'s hydrogen-envelope branch already puts the surface on a dry
adiabat and checks the boiling curve at envelope pressure. Add `hycean` to the
regimes when an ocean sits under an H2 envelope, below the boiling curve at that
pressure and shallower than the 5.2 limit. Deeper ones are ice-sealed
([atmospheres-retention-and-classes.md](atmospheres-retention-and-classes.md)
4.3).

## 6. Corrections to the source documents

| # | Source and section | Statement | Correction |
|---|---|---|---|
| B1 | Rad "Host Star Spectral Regimes" and its table | `L_XUV/L_bol` about 1e-3 saturated for M, K, G | That is the X-ray band. Total XUV is 3.3e-3 (G) to 7.3e-3 (0.1 Msun); F stars saturate lower (down to 10^-4.3) [C] |
| B2 | Rad table | M-dwarf `t_sat` 1.0 to 2.5 Gyr; G-dwarf flares over 1e34 erg 1e-4 to 1e-2 per year | `t_sat` fine for 0.3 to 0.5 Msun, 3 to 8 Gyr for late M; flares 1.7e-4 to 1.2e-3 per year for slow rotators, 1e-2 only for young fast ones |
| B3 | Rad "Atmospheric / Magnetic Regime" table; Master radiation table | "Earth-equivalent about 2.4 mSv/yr" as surface dose; quiescent GCR 2.4 mSv/yr | 2.4 is the total natural background; cosmic is 0.39 mSv/yr |
| B4 | Rad "Planetary Shielding Dynamics" | A dipole reduces dose by two or three orders of magnitude | True for solar particles only. GCR: 0.6 at most (4.2) |
| B5 | Rad table, "Rarefied Unmagnetized (0.01 bar)" | More than 10 to 100 Gy per flare | Needs an extreme flare. For the Sun a 0.01 bar atmosphere gives about 2 mSv/yr of SEP plus about 265 mSv/yr of GCR [C] |
| B6 | Rad "K-Dwarfs" | K-dwarf habitable-zone planets avoid tidal locking | Holds for 0.75 to 0.85 Msun only; locked by 10 Gyr below about 0.7 Msun (3.3) |
| B7 | Rad "M-Dwarfs" | PMS contractions of 1.0 to 2.5 Gyr with saturated XUV | Contraction and saturation are different clocks; see A10 in [atmospheres-retention-and-classes.md](atmospheres-retention-and-classes.md) |
| B8 | PHI-4 Radiation; Master `L_rad`; habitability-index.md (old) | `L_rad` threshold 10 Gy/yr disagrees with the Mars row (`L_rad` 0.65) | A unit slip (mSv against Gy). Mars is 0.077 Gy/yr, 130x below 10, so `L_rad(Mars)` is about 1 |
| B9 | PHI-4 Radiation tiers | Blue < 50 mSv/yr, Green 50 to 100, Yellow "subsurface", Red > 10 Gy/yr | Gap from 100 mSv/yr to 10 Gy/yr; filled in 4.6 |
| B10 | PHI-4 Radiation; Master `M_rad` | Blue 50 mSv/yr (occupational) | ICRP: 20 mSv/yr (five-year average), 50 in any one year. Default Blue 20 |
| B11 | Chem "Hyperacidic Sulfate Oceans" and table; PHI-4 Chemistry Red; Master | pH ranges differ (1.0 to 4.0; pH 0 to 1) | Use 1.0 to 4.5 for the class; pH 0 to 1 is the tail. Soda: Chem 9.0 to 11.5, Master table 9.5 to 11.5; use 9.0 to 11.5 |
| B12 | Chem "High-Pressure Ice-Sealed Oceanworlds"; Master; Math 2.1 | Ice VI/VII seal above 1 to 1.2 GPa | Depends on temperature: 0.2 GPa near 251 K, 0.63 GPa at 273 K, 2.2 GPa at 355 K; 73 to 204 km of ocean at Earth gravity (5.2) |
| B13 | Chem "Thermodynamic Extremes"; PHI-4 Temperature Yellow | Moist-greenhouse surface 50 to 90 C | The moist limit is a flux threshold (`S_eff` 1.01 for the Sun); the surface temperature at that flux depends on the atmosphere [R]. Zone edges are in [atmospheres-retention-and-classes.md](atmospheres-retention-and-classes.md) 4.1. The Chem passage also has a garbled LaTeX run in the silica sentence |
| B14 | Math 2.4 | Rigidity cutoff and Chapman-Ferraro standoff named without formulas | Formulas in 3.4 and 4.2 |
| B15 | Math 1.1 Category A, rotation | Rotation from a log-normal unless `t_sys >= t_sync` | Agrees with 3.3; the project draws 10 to 1,400 h uniformly and locks only moons |
| B16 | rogue-planet-surface.md "Simplifications" | Ignores high-pressure ice | Must be fixed once ocean depth feeds the chemistry class (5.3) |

## 7. Evidence notes

Items tagged [S] were confirmed in a search summary only. [R] claims to verify
when paper access is allowed:

1. Saturated `L_X/L_bol` = 10^-3.13 (one search summary only).
2. The Olson and Christensen 2006 equation, gamma = 0.2 and the `F`
   definition; the `Ro_l` = 0.12 threshold and Earth's 0.09.
3. The Sanz-Forcada/Chadney `L_EUV` relation (4.80, 0.860) and Ribas's
   4.6 erg/cm2/s; Wright 2018's slope for fully convective stars; Wood 2005
   wind scaling and the cap values.
4. Young solar-type flare rates (all non-anchor cells of the 2.4 table) and the
   SEP scaling (20 mSv/yr, the square-root exponent, `R0` = 0.2 GV from the
   1956 event, the 5.4 Sv/day Europa figure).
5. GCR dose at 100 g/cm2, the 164 g/cm2 attenuation length and the heliosphere
   multiplier (2.5x).
6. Peale's tidal formula and k2/Q values; ice-shell convective onset; Earth
   ocean pH and salinity; the dynamo lifetime rules; ICRP limits.

## 8. Sources

Search results used (the full papers could not be opened):

- Activity: https://arxiv.org/pdf/2010.12922, https://arxiv.org/pdf/1502.07401,
  https://ar5iv.labs.arxiv.org/html/1111.0031,
  https://ar5iv.arxiv.org/html/0812.1223, https://arxiv.org/abs/astro-ph/0412253,
  https://www.arxiv.org/abs/1504.04546, https://arxiv.org/pdf/1412.3380
- Flares: https://arxiv.org/pdf/1904.00142, https://ar5iv.labs.arxiv.org/html/2011.02117,
  https://arxiv.org/pdf/1912.11572, https://arxiv.org/pdf/1804.02001,
  https://arxiv.org/pdf/1904.06875, https://arxiv.org/pdf/1907.12580,
  https://www.aanda.org/10.1051/0004-6361/202142710, https://arxiv.org/pdf/1901.00890
- Ozone and flares: https://arxiv.org/abs/1006.0022 (Segura),
  https://arxiv.org/pdf/1711.08484
- Dynamo and magnetopause: https://www.maths.gla.ac.uk/~rs/res/B/DynamoScaling/2006OlsonChristensen.pdf
  (search result only), https://arxiv.org/pdf/1302.7140,
  https://ar5iv.labs.arxiv.org/html/1101.0691, https://arxiv.org/pdf/1304.2909,
  https://arxiv.org/pdf/0902.0952, https://arxiv.org/pdf/1509.00735
- Dose: https://elib.dlr.de/106090 (MSL RAD), https://arxiv.org/pdf/2001.11028
  (Chang'E-4 LND),
  https://www.env.go.jp/en/chemi/rhm/basic-info/2018/pdf/basic-1st-02-05-slides.pdf
  (UNSCEAR 2008), https://unis.unvienna.org/unis/en/pressrels/2026/unisous453.html,
  https://icrpaedia.org/LD_50/30,
  https://publish.nrc.gov/reading-rm/basic-ref/glossary/lethal-dose-ld.html,
  https://intechopen.com/chapters/48814,
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5173286/
- Supernovae and GRBs: https://ar5iv.arxiv.org/html/astro-ph/0211361 (Gehrels),
  https://ar5iv.labs.arxiv.org/html/astro-ph/0505472 (Thomas),
  https://arxiv.org/pdf/astro-ph/0411284.pdf, https://arxiv.org/pdf/1107.1286
  (Gowanlock), https://arxiv.org/pdf/2009.13539.pdf (Spinelli)
- Hydrosphere: https://arxiv.org/pdf/1801.00748 (Kite and Ford),
  https://www.nature.com/articles/nphys3822 (Noack),
  https://arxiv.org/abs/1401.0720 (Cowan and Abbot),
  https://arxiv.org/pdf/2108.10888 (Hycean),
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12827049/ (Juno, Europa),
  https://arxiv.org/pdf/2503.01967 (Enceladus), https://arxiv.org/pdf/1803.01511
  (Ganymede), https://arxiv.org/pdf/1712.03614, https://arxiv.org/pdf/1903.12111,
  https://arxiv.org/pdf/1907.09598 and https://water.lsbu.ac.uk/water/water_phase_diagram.html
  (phase diagram)
