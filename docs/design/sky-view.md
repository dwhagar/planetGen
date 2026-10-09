# The view from a planet: sky rendering and constellations

How to compute and draw the sky a planet, moon or facility sees: which stars are visible and how bright, what colour, where each one appears (light-time, parallax, the planet's own orientation), how the galaxy's unresolved glow is added, how constellations are chosen and named, and how the picture is produced and cached. It carries the math Boss's VIEW.1 asked for before any VIEW.2 code is written, checked against a real star catalogue (HYG v4.1) and against the repo's code. The other galaxies in that sky are in [multiple-galaxies.md](multiple-galaxies.md).

Informs: VIEW.1, VIEW.2, VIEW.3, VIEW.4 (VIEW.5, GEN.6 and GEN.125 are done; GEN.104 is a prerequisite for horizon views)

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [S] seen in a search result or a file fetched from raw.githubusercontent.com, [C] computed in this research (numpy and scipy against HYG v4.1, the repo's own functions, the CIE tables in `colour-science`), [R] recalled and not confirmed. The research environment could read search-result text and raw GitHub files but not the papers, so every [R] item is listed in the evidence notes at the end. The scripts were run from a scratch folder and are not in the repo; the function bodies below are the reference and reproduce the quoted numbers.

## Summary

- **The naked-eye sky is small and local.** HYG v4.1 has about 8,920 stars at V <= 6.5 [C], against the Yale catalogue's 9,096 [S]; 88% lie within 300 pc and 99% within 750 pc [C]. A planet's sky needs about 9,000 resolved stars plus a statistical glow, not millions.
- **The existing bright-star tables cover about 8% of it.** The pre-placed `bright_stars` (floor 1,000 L_sun) and the 100 ly backfill leave most naked-eye stars unmade; 75% needs L >= 30 L_sun out to about 110 pc. The plan adds sky tiers on the existing deterministic band backfill (section 3).
- **Colour is cheap.** A blackbody integrated through the CIE curves matches `_kelvin_to_hex` within about 8 units per channel [C]; add dust reddening and desaturate faint stars (section 4).
- **Light-time barely moves stars** (median 15 arcsec [C]); what matters is the star's age and luminosity at the retarded time. **Parallax moves them a lot**: 20 ly shifts the median naked-eye star by 2.2 degrees [C], so every system needs its own sky and constellations.
- **Constellations** need a partition (seeds 20 degrees apart, one MST per region, 12 stars each), not one global MST: 56 figures over the 600 brightest stars, median segment 5.8 degrees against the real 4.8 (section 6).
- **Rendering** needs no new install: numpy and `zlib` write the PNG, and Pillow is already present as scikit-image's dependency (section 7).
- **Prerequisite:** GEN.104 (spin axis, tilt, rotation phase) for horizon views and time of day. Whole-sky charts do not need it.

## Decisions already taken

- Boss (VIEW.1, 2026-10-01, TODO items 83 to 85): "view-from-planet will have to do calculations on colors and A LOT Of stuff, so make special note of that, it will need a full research pass." VIEW.2 and VIEW.3 wait for a research session with him; VIEW.4 does not.
- Boss (VIEW.2): the view places each star "at that light years back in time", that is where it was when the light now arriving left it. Done as VIEW.5 (PR #719), with the galactic orbits of GEN.6 (7.37.0, PR #157).
- Boss (VIEW.3): the view is drawn to a PNG with constellations, stored so a planet keeps the same constellations each time.
- Plan of 2026-10-07 and Boss, 2026-10-08 03:57Z and 04:00Z ([object-ids.md](object-ids.md)): constellation names come from the codec under the naming key, not from sliced word lists. `names/naming_key.py` already lists `constellation` in `KINDS`.
- Boss, 2026-09-30 ([navigation-frames.md](navigation-frames.md)): east is North x Up. The horizon frame below uses the same construction.

## 1. What exists

| Item | Status | Used here for |
|---|---|---|
| VIEW.5 `physics/light_travel.py` | done, PR #719 | `apparent_position(position_at, observer, now_s)` iterates the retarded time; `apparent_position_of(body, observer)` for a `SpatialPosition3D`. One body at a time; the sky needs the vectorised version in section 2 |
| GEN.6 galactic motion | done, 7.37.0, PR #157 | stars and systems move along galactic orbits |
| GEN.125 facility velocity | done, PR #782 | a stand-alone facility is an observer with a position and a velocity |
| GEN.98 backfill from the farthest boundary outward | done, PR #795 | the tier machinery of section 3 |
| GEN.104 spin vector and axial tilt | open | horizon, time of day, seasons |

Data the sky reads:

- **Positions:** systems sit at `position_*_mpc` inside their sector (sector-local frame, [galaxy-coordinate-system.md](galaxy-coordinate-system.md) section 2); `bright_stars` has integer `position_*_mpc`. The 1 mpc quantisation is 206 AU, 0.02 degrees at 1.3 pc and irrelevant beyond.
- **Velocities:** `star_systems.velocity_*_kms` (v61, GEN.121), planets and moons (v60), facilities (v64). Pre-placed `bright_stars` rows carry none until their system is built; derive one from the rotation curve (`galaxy.system_position.galactic_velocity_ms`) if needed.
- **Stars store** `luminosity_w`, `temperature_k`, `radius_km`, `age_gy`, `lifespan_gy` and (bright stars) `phase_end_age_gy`. No colour index and no magnitude: compute both on the fly; no schema change.
- **Dust:** `nebulae.extinction_av`, `radius_ly`, `nebula_shape_balls`; `supernova_remnants.age_years`, `radius_ly`.
- **Density:** `galaxy/density.py::relative_density` and `model_terms(shape)['solar_angle_rad']`. The model places "the Sun" at R0 = 8.2 kpc on the inter-arm minimum (azimuth 5.8575 rad for the default shape [C]); the bar's near end leads it by 27 degrees ([galaxy-disk-density.md](galaxy-disk-density.md)).
- **Colour code today:** `starmap.py::star_color` and `_COLOR_NAME_RGB` use a named colour per spectral letter (hand-picked, the file says "not colorimetric"); `_kelvin_to_hex` is the Tanner Helland fit, used as a fallback.
- **Rotation sense:** the project's galaxy rotates counterclockwise seen from +Z, the real one clockwise seen from its north pole ([galactic-potential.md](galactic-potential.md) section 5). Any real-sky data (a catalogue star list, the Local Group) needs a mirror or the flipped embedding in [multiple-galaxies.md](multiple-galaxies.md) section 2.

## 2. The math pipeline

Units: parsec, year, solar luminosity, kelvin.

**2.1 Observer and vectors.** `P_obs` is the planet's system position in the galactic frame (sector centre plus system offset, rotated by `sector_orientation`) plus the planet's orbital offset (4.85e-6 pc per AU). The orbital offset is irrelevant for most stars and the whole point for the planet's own star and companions. For each star `v = P_star - P_obs`, `d = |v|`, `u = v/d`.

**2.2 Light-time.** Replace `P_star` by `P_star - V_star * tau` and iterate twice with `tau[yr] = d[pc] * 3.2616`. Each pass shrinks the error by about v/c = 1e-4, so two passes are exact to the float. Straight-line motion is exact enough out to kiloparsecs.

```python
import numpy as np
LY_PER_PC = 3.261563777
MBOL_SUN = 4.74

def sky_vectors(star_pc, observer_pc, star_vel_pc_yr=None):
    """Unit vectors and distances to stars, with first-order light-time (pc, pc/yr)."""
    v = star_pc - observer_pc
    d = np.linalg.norm(v, axis=1)
    if star_vel_pc_yr is not None:                       # seen where it was d/c years ago
        for _ in range(2):
            v = star_pc - star_vel_pc_yr * (d * LY_PER_PC)[:, None] - observer_pc
            d = np.linalg.norm(v, axis=1)
    return v / d[:, None], d
```

**2.3 Light-time in age and luminosity.** The star is seen at age `a' = age - tau`. A star whose lifespan ended less than `tau` ago is still seen alive; for a supernova remnant of age `A` at distance `d`, the explosion is visible only if `A > d/c`, otherwise the sky shows the progenitor. The stars affected have lifetimes under about 10 Myr (a light-travel of 3 kpc); the naked-eye ones sit at 1 to 3 kpc (an O star with M_V = -6 reaches V = 6.5 at 3.2 kpc). A correctness item rather than a visual one.

**2.4 Magnitude.** `M_bol = 4.74 - 2.5 log10(L/L_sun)`, `M_V = M_bol - BC_V(T)`, `m = M_V - 5 + 5 log10(d_pc) + A_V`. BC_V comes from integrating a blackbody through a V-like response, calibrated to the Sun's -0.07 [C]:

| T (K) | 2,500 | 3,000 | 3,500 | 4,000 | 5,000 | 5,772 | 7,000 | 10,000 | 20,000 | 30,000 | 45,000 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BC_V | -2.72 | -1.69 | -1.05 | -0.63 | -0.20 | -0.07 | -0.04 | -0.32 | -1.64 | -2.68 | -3.82 |

These agree with remembered main-sequence values to about 0.3 mag [R]; the hot end is weakest (the real BC at 45,000 K is about -4.4 [R]). A 0.3 mag error changes the count at V <= 6.5 by about 15%, so replace the table with a published relation (Flower 1996 or Torres 2010) before anything depends on the count.

```python
_T = np.log([2500, 3000, 3500, 4000, 5000, 5772, 7000, 10000, 20000, 30000, 45000])
_BC = np.array([-2.72, -1.69, -1.05, -0.63, -0.20, -0.07, -0.04, -0.32, -1.64, -2.68, -3.82])
def apparent_mag(lum_sol, teff_k, d_pc, a_v=0.0):
    bc = np.interp(np.log(teff_k), _T, _BC)
    return MBOL_SUN - 2.5 * np.log10(lum_sol) - bc + 5 * np.log10(d_pc / 10.0) + a_v
```

Pitfall: HYG's `lum` column is V-band luminosity, not bolometric. Feeding it to the bolometric formula double counts BC and gave 6,651 stars at V <= 6.5 instead of 8,639 [C]. The project's `luminosity_w` is bolometric, so the formula applies as written.

**2.5 Dust.** `A_V` is the dust density integrated along the line of sight, in three tiers: (a) the chord through each `nebulae` and remnant ball, scaled by `extinction_av` over the chord-to-radius ratio; (b) a smooth galactic law `a0 exp(-R/R_d) exp(-|z|/h_d)`; (c) nothing, for early versions. Averages are 0.7 to 1.0 mag/kpc near the Sun and about 1.8 mag/kpc in the plane over a few kpc [S, an encyclopedia mirror of Wikipedia "Extinction (astronomy)"; primary source not confirmed]. Default for (b): `a0` set so the local plane gives 1.0 mag/kpc, `R_d = 2.26 kpc`, `h_d = 134 pc` [R: Drimmel and Spergel 2001]. Reddening: `E(B-V) = A_V/3.1`, with the spectrum multiplied by `10^(-0.4 A_V (a(x) + b(x)/R_V))` using the Cardelli, Clayton and Mathis (1989) optical polynomials. The recalled coefficients give `A_B/A_V = 1.3245` against the defining 1.3226 for R_V = 3.1 [C]. A 5,772 K star at E(B-V) = 0.3 moves from #fff1ea to #ffe2c3; a 10,000 K star at E(B-V) = 1.0 from #cdd9ff to #ffdbad [C].

**2.6 Planet orientation (needs GEN.104).** Inputs: the spin axis `s` (unit vector, galactic frame) and the sidereal rotation rate. The kinetics document ([Observational Kinetics for Rotational Vectors.md](Observational%20Kinetics%20for%20Rotational%20Vectors.md), "IAU Cartographic Coordinate Systems") gives the standard form: pole `(alpha0, delta0)`, prime meridian `W(t) = W0 + Wdot d` with `Wdot = 360 / P_rot`, and `omega = (cos d0 cos a0, cos d0 sin a0, sin d0)`. [orbital-updates.md](orbital-updates.md) section 6 and [orbital-solvers-and-integrators.md](orbital-solvers-and-integrators.md) section 9 take the pole in the galactic frame.

Two checks against the kinetics document:

- Its expanded 3x3 matrix does not equal the product it states, `Rz(90 + a0) Rx(90 - d0) Rz(W)`. Only the third column (the pole) agrees; the first two columns of the expansion equal the product's with W shifted by 90 degrees [C, numerical check at a0 = 30, d0 = 40, W = 50 degrees]. Use the product. At W = 0 it puts the body's x axis at `(-sin a0, cos a0, 0)`, the ascending node of the body's equator on the reference plane, which is the prime-meridian reference direction below.
- The document's cascade makes the sky's geometry follow from three outputs of GEN.104: the axis `s`, a sidereal period, and `W0` at an epoch. A tidally locked body (`P_rot = P_orb`, obliquity 0) has a fixed sun and a parent that hangs at one altitude and azimuth; the stars wheel once per orbit. An obliquity above 90 degrees (the planet-collision distribution allows it) reverses the sense: the sun rises in the west. The formulas below need no special case for either.

The generator today sets `rotation_period_hours` from a range and, for a tidally locked moon, to the orbital period (`physics/planets.py`), which is a sidereal period (the solar day of a locked body is infinite). The column therefore reads as sidereal. It is documented as "a static descriptive stat; no rotational phase is tracked" ([database-schema.md](../database-schema.md)), so GEN.104 must also store `W0` or derive it from the seed, or time of day has no meaning. Only moons are tidally locked today; close-in planets are not.

Body frame: `z = s`, reference direction `r = (-sin a0, cos a0, 0)` (or any seeded direction), `x = normalize(r - (r.s) s)`, `y = z x x`. An observer at latitude `phi` and longitude `lam` with rotation angle `W = W0 + 360 t/P` has the zenith `U = cos(phi)(cos(lam+W) x + sin(lam+W) y) + sin(phi) z`. North is `N = normalize(s - (s.U) U)` (at a pole, any horizontal direction), east `E = N x U`. Altitude is `asin(v.U)`, azimuth `atan2(v.E, v.N)`. This is the construction of [navigation-frames.md](navigation-frames.md) with the spin pole in place of "toward the centre of the frame", so the planet's N and the navigation N are different things; the sky page must say which it uses.

```python
def zenith_inertial(spin_axis, ref_dir, lat_deg, lon_deg, rot_angle_deg):
    z = spin_axis / np.linalg.norm(spin_axis)
    x = ref_dir - (ref_dir @ z) * z; x /= np.linalg.norm(x); y = np.cross(z, x)
    lon = np.radians(lon_deg + rot_angle_deg); lat = np.radians(lat_deg)
    return np.cos(lat) * (np.cos(lon) * x + np.sin(lon) * y) + np.sin(lat) * z

def horizon_basis(spin_axis, zenith):
    up = zenith / np.linalg.norm(zenith)
    n = spin_axis - (spin_axis @ up) * up
    n = n / np.linalg.norm(n) if np.linalg.norm(n) > 1e-9 else np.cross(up, [1.0, 0, 0])
    return n, np.cross(n, up), up                        # north, east (N x U), up
```

The planet's own star comes from the orbit (`physics/body_positions.py`). Day and night follow from `U` and the star direction: sun altitude above 0 is day, -6 to -18 degrees is twilight, and the limiting magnitude deepens only after astronomical twilight ends (-18). Axis precession is a ten-thousand-year effect; ignore it.

**2.7 Atmosphere and limiting magnitude.** The naked-eye limit is about 6.5 [S: the Yale catalogue lists stars to 6.5, 9,096 of them, about half visible at once]. Near the horizon `m' = m + k X(alt)`, airmass `X ~ 1/(sin(alt) + 0.025 exp(-11 sin(alt)))`, `k` about 0.2 mag per airmass for a clear Earth sky [R]; scale `k` with surface pressure. A world with a thick cloud deck has no sky and returns an explicit "no view".

**2.8 Aberration** of the observer at 220 km/s is 0.042 degrees [C]; ignore.

**2.9 Planets, moons and companions.** The planet's own star and any companion use the same formula with `d` in AU (m = -26.7 at 1 AU for a Sun-like star). Other bodies in the system shine by reflection: star flux at the body times `albedo * (R/d)^2 * phase(alpha)` with a Lambert-type phase function, converted with the same `-2.5 log10` rule, drawn as wandering points. Out of scope for the first stage.

**2.10 Extended objects.** Nebulae: an ellipse sprite with total magnitude and angular radius `atan(radius_ly/d)`. Other galaxies: [multiple-galaxies.md](multiple-galaxies.md). Visibility is a surface-brightness test: an object is naked-eye only if its mean surface brightness is below about 23 mag/arcsec^2 [R] (M31 22.2, M33 23.0, the LMC 22.6, the SMC 23.3 by calculation from catalogue values [C]).

## 3. Which stars to generate: sky tiers

### 3.1 Resolved set and glow

Resolved stars come from the real generator; everything fainter is integrated light. Never invent resolved stars for unfilled space, because a later sector fill would not reproduce them.

1. **Resolved set** (V <= m_lim, default 6.5): `bright_stars`, generated sectors, and sky-tier backfills (3.3).
2. **Unresolved glow:** for each sight line `I = (1/4pi) integral j(s) f_faint(s) 10^(-0.4 A_V(s)) ds`, with `j(s) = j0 * relative_density(P_obs + s u)`, `j0` the local V-band luminosity density, and `f_faint(s)` the share of the luminosity function dimmer than the star that would reach `m_lim` at distance `s` (so resolved stars are not counted twice). Measured `j0 = 0.044 L_V,sun/pc^3` from the 3,081 HYG stars within 25 pc (incomplete for M dwarfs, which carry little light; the literature is about 0.05 [R]). In that sample stars dimmer than 1 L_sun carry 7% of the local V light, dimmer than 10 L_sun 35%, dimmer than 100 L_sun 81% [C]. Surface brightness is `SB = -0.19 - 2.5 log10(4 pi I Omega_arcsec2)`; omitting the 4 pi made the first run's Milky Way 3 mag too faint, so test this against known values.
3. **Pre-rendered band:** ray-march the project's own density from the observer. The numpy mirror of `relative_density` agrees with the Python function to 1.9e-15 relative on 2,000 random points [C]. A 180x90 map with 900 steps takes 13 s, 72x36 with 500 steps about 2 s. For the default shape at the Sun with 1.0 mag/kpc dust: 21.6 mag/arcsec^2 at (l = 0, b = +2.5), 21.8 at the anticentre, 23.0 at the pole [C], the right range for the Milky Way's naked-eye band (about 21 to 23 [R]). With no dust the bulge direction reaches 15.8, so the dust model decides how the galactic centre looks: tune it. Moving the observer 20 ly changes the map by at most 0.004 mag, so cache one band per sector or per 50 pc cell, not per planet. From R = 14 kpc the in-plane range is 18.5 to 23.3; from R = 500 pc the sky is bulge-lit everywhere (pole 20.6) [C].
4. Demo images made from real HYG stars plus the band (Hammer-Aitoff 1600x800 in 4 s, stereographic horizon 1000x1000 in 1.5 s) looked plausible: bulge, dust rift, colours. Flaws to fix: 1-degree band pixelisation (use bilinear lookup and a finer map) and 8-bit contouring in the glow (add dither).

### 3.2 Cost model

Naked-eye stars by absolute magnitude (HYG v4.1, V <= 6.5, known distance, 8,714 stars) [C]:

| M_V bin | stars | median d (pc) | max d (pc) |
|---|---:|---:|---:|
| -5 and brighter | 328 | about 570 | 980 |
| -3 to -1 | 1,561 | 301 | 763 |
| -1 to 1 | 4,143 | 135 | 310 |
| 1 to 3 | 2,025 | 66 | 125 |
| 3 to 5 | 576 | 28 | 48 |
| 5 and fainter | 81 | about 10 | 19 |

Share of the naked-eye sky above a luminosity floor (M_V from L with BC = -0.3): L >= 1,000 L_sun 8%, >= 250 27%, >= 100 42%, >= 30 75%, >= 10 89% [C]. Only 602 naked-eye stars lie within 100 ly.

The deterministic machinery exists. `generation/bright_stars.py` splits luminosity into eight fixed bands per decade, each cell's band drawn from its own stream (`canonical_bands`), so stars between two levels never depend on the steps taken to get there. `sector_stats.bright_level_sol` records how deep a cell has been drawn, and `tuning.BRIGHT_STAR_BACKFILL_TIERS` is `((10, 100), (25, 250), (50, 500), (100, 750))` as `(out_to_ly, min_luminosity_sol)`. A sky needs a second tier table whose floors follow visibility. Cells are 4 pc cubes; the naked-eye distance for a floor L is `10^((6.5 - M_V + 5)/5)` pc, ignoring dust:

| L floor (L_sun) | naked-eye radius | cells within | radius at V <= 7.5 | cells |
|---:|---:|---:|---:|---:|
| 1,000 | 619 pc (2,020 ly) | 1.6e7 | 982 pc | 6.2e7 |
| 300 | 339 pc | 2.6e6 | 538 pc | 1.0e7 |
| 100 | 196 pc (639 ly) | 4.9e5 | 310 pc | 2.0e6 |
| 30 | 107 pc (350 ly) | 8.1e4 | 170 pc | 3.2e5 |
| 10 | 62 pc (202 ly) | 1.6e4 | 98 pc | 6.2e4 |
| 3 | 34 pc | 2.6e3 | 54 pc | 1.0e4 |
| 1 | 20 pc | 490 | 31 pc | 2.0e3 |

The project's whole backfill (100 ly) is 1,886 cells. The 1,000 L_sun level is already galaxy-wide in `bright_stars` (26.9 million stars in the 2026-10-08 timing run, `bright-star-timing/report.md`); 100 to 300 L_sun is the expensive middle, 0.5 to 2.6 million cells per sky, where a cell-by-cell Poisson draw is mostly empty. The backfill's cost per cell was measured later: about 0.06 ms per cell for the draws, about 2 ms per sector once warm and about 5 s for a cold start ([sampling-backfill-and-resume.md](sampling-backfill-and-resume.md)); the earlier timing run measured the sector fill, 46.6 ms per bright star, not `backfill_cells`. Two ways to cut it: (a) share the work, since `sector_stats` remembers each cell's level and neighbouring planets reuse it; (b) draw the middle bands ring by ring in angle bins as `scatter_layer` does for the whole galaxy, keeping the same per-band streams. Dust helps: at 1 mag/kpc the 1,000 L_sun radius drops from 619 to about 480 pc.

### 3.3 Rule

A sky is built from (i) `bright_stars` and generated stars within the sphere, (ii) a backfill to the sky tiers for the planet's system, (iii) the band for the rest. A sky is complete for a limiting magnitude `m_lim` only once its tiers have run. Store the tier level reached (`sector_stats.bright_level_sol` already does) and return "sky incomplete" or run the backfill.

## 4. Colour

Planck spectrum at `T_eff` into XYZ with the Wyman, Sloan and Shirley (2013, J. Computer Graphics Techniques 2(2)) multi-lobe Gaussian fit of the CIE 1931 curves, into linear sRGB with the D65 matrix, scaled so the maximum channel is 1. The recalled Gaussian coefficients match the tabulated CIE 1931 2-degree observer with a maximum absolute error of 0.015, 0.007 and 0.024 against peaks 1.06, 1.00 and 1.78 [C].

```python
def _g(l, mu, s1, s2): return np.exp(-0.5 * ((l - mu) / np.where(l < mu, s1, s2)) ** 2)
def blackbody_rgb(T, lam=np.arange(380.0, 781.0)):
    l = lam * 1e-9
    s = 1 / (l ** 5 * np.expm1(1.438777e-2 / (l * T)))
    x = 1.056*_g(lam,599.8,37.9,31.0) + 0.362*_g(lam,442.0,16.0,26.7) - 0.065*_g(lam,501.1,20.4,26.2)
    y = 0.821*_g(lam,568.8,46.9,40.5) + 0.286*_g(lam,530.9,16.3,31.1)
    z = 1.217*_g(lam,437.0,11.8,36.0) + 0.681*_g(lam,459.0,26.0,13.8)
    XYZ = np.array([(s * c).sum() for c in (x, y, z)])
    M = np.array([[3.2404542,-1.5371385,-0.4985314],[-0.969266,1.8760108,0.041556],[0.0556434,-0.2040259,1.0572252]])
    rgb = np.clip(M @ XYZ, 0, None); return rgb / rgb.max()
_TG = np.geomspace(2200, 45000, 96); _RGB = np.array([blackbody_rgb(t) for t in _TG])
def star_rgb(teff_k):
    lt = np.log(np.asarray(teff_k)); return np.stack([np.interp(lt, np.log(_TG), _RGB[:, k]) for k in range(3)], -1)
def srgb8(lin):                                          # linear -> 8-bit sRGB; do this after summing light, not before
    lin = np.clip(lin, 0, 1); return np.round(255 * np.where(lin <= 0.0031308, 12.92 * lin, 1.055 * lin ** (1 / 2.4) - 0.055)).astype(int)
```

`star_rgb` returns linear RGB; the hex values in the table below are `srgb8` of it.

| T (K) | computed sRGB | project `_kelvin_to_hex` |
|---:|---|---|
| 3,000 | #ffb96e | #ffb16e |
| 4,000 | #ffd4a6 | #ffcea6 |
| 5,772 | #fff1ea | #fff2e6 |
| 7,500 | #ebecff | #e6ebff |
| 10,000 | #cdd9ff | #cadaff |
| 30,000 | #a2bbff | #9fbeff |

Rules:

- The hand-picked anchors in `_COLOR_NAME_RGB` are a map style, not sky colour; do not reuse them for the sky. `_kelvin_to_hex` may stay for the maps.
- **Perception.** Faint stars look white to the eye (rods) and only the brightest, such as Betelgeuse and Antares, show colour [S, via search: CSIRO, Gresham College, Astronomy magazine]. Blend toward white by magnitude: `sat = clip(0.35 + 0.65 (5.5 - m)/4, 0.35, 1)`, `rgb' = 1 - (1 - rgb) sat`. The constant is a starting guess; tune by eye. Offer a "boost colour" option for the stronger map look.
- **B-V to T_eff,** needed only when importing a real catalogue: Ballesteros (2012, EPL 97, 34008), `T = 4600 (1/(0.92 (B-V) + 1.7) + 1/(0.92 (B-V) + 0.62))`; the Sun's B-V 0.65 gives 5,778 K [C]. The generator has `temperature_k` directly.
- Dust reddening is section 2.5. Atmosphere colour, the planet's own star's glare and twilight tint belong to the horizon stage.

## 5. Projections

Formulas for a unit vector at angle `theta` from the view centre; `r` is the image radius. Forward-then-inverse round trips agree to 2e-6 degrees on 200,000 random directions for each, and Hammer is equal-area (counts per equal-area ring flat within 1%) [C].

| Projection | Formula | Property | Use |
|---|---|---|---|
| Stereographic | `r = 2 tan(theta/2)` | conformal; linear scale 1.33 at 60 degrees, 2.0 at 90, 4.0 at 120 [C] | the horizon dome up to 90 to 120 degrees; constellations keep their shape |
| Azimuthal equidistant | `r = theta` | true distance from the centre | whole sky as a disc |
| Equisolid fisheye | `r = 2 sin(theta/2)` | equal-area | camera-like fisheye |
| Orthographic | `r = sin(theta)` | one hemisphere | the globe look |
| Hammer-Aitoff | `x = 2 sqrt2 cos(lat) sin(lon/2) / d`, `y = sqrt2 sin(lat) / d`, `d = sqrt(1 + cos(lat) cos(lon/2))`; inverse `z = sqrt(1 - (x/4)^2 - (y/2)^2)`, `lon = 2 atan2(z x, 2(2 z^2 - 1))`, `lat = asin(z y)` | equal-area, whole sky in an ellipse | the all-sky chart and the band |
| Equirectangular | `x = lon`, `y = lat` | simple, distorts the poles | texture storage |

The inverse Hammer `atan2` needs the factor 2 in its second argument; without it the round trip is off by 39 degrees.

Recommendation: whole sky as Hammer-Aitoff centred on the galactic centre direction (l on the horizontal axis, increasing to the left as in catalogues); horizon view as stereographic centred on the zenith with the horizon circle at 90 degrees, field of view up to 150 degrees if wanted.

```python
def stereographic(vec, basis):
    n, e, u = basis; vx, vy, vz = vec @ e, vec @ n, vec @ u
    th = np.arccos(np.clip(vz, -1, 1)); r = 2 * np.tan(th / 2); ph = np.arctan2(vy, vx)
    return np.stack([r * np.cos(ph), r * np.sin(ph)], 1), vz > 0
```

## 6. Constellations

### 6.1 Measured

Real IAU figures (89 polylines from d3-celestial's `constellations.lines.json`, Serpens counted twice) [C]: median 7 stars and 7 segments, median segment 4.8 degrees (p10 1.9, p90 10.1, max 25.6), 1.00 segments per star, about 700 stars down to V = 4.4.

| Method over the N brightest stars | N | figures (3+ stars) | median segment | note |
|---|---:|---:|---:|---|
| Delaunay (spherical, via convex hull) | 100 | 1 | 20 deg | 2.9 edges per star, far too dense |
| MST | 100 | 1 | 10.2 deg | one tree |
| MST, edges over 20 deg cut | 100 | 7 | 9.0 deg | |
| RNG (relative neighbourhood), cut 25 deg | 100 | 5 | 9.6 deg | 1.0 edge per star |
| Seeded regions, 30 deg, cap 12 | 100 | 15 | 8.9 deg | |
| MST cut 8 deg | 600 | 43 | 4.2 deg | many 2-star scraps |
| MST cut 12 deg | 600 | 3 | 4.8 deg | percolates into a few giants |
| RNG cut 12 deg | 900 | 1 | 4.7 deg | percolates |
| **Seeded regions, 20 deg, cap 12** | 600 | **56** (508 stars) | **5.8 deg** | closest to the IAU shape |

A global MST or RNG percolates once the star set is dense enough to match real segment lengths, and cutting by length shatters the figures at the threshold that stops the percolation. The brief's "brightest ~100 stars" gives segments of about 10 degrees and is too sparse: use the 500 to 900 brightest.

### 6.2 Algorithm

1. Take the `N` brightest visible stars, sorted by `(m, star uid)` so ties are deterministic. `N` about 600 (V <= 4.1 from Earth); scale `N` to the sky so figure stars are about one per 70 square degrees.
2. Seeds: walk the list brightest first; a star becomes a seed if it is at least 20 degrees from every earlier seed.
3. Every star joins its nearest seed. This partition gives every star exactly one constellation, so Bayer lettering is well defined (the IAU constellations also partition the sphere).
4. In each region keep the 12 brightest stars and join them with the MST over their spherical Delaunay edges (a subset of Delaunay, so edges never cross). Optionally add one or two short RNG edges to close a loop in a bare chain.
5. Drop regions with fewer than 3 stars.

```python
from scipy.spatial import ConvexHull
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import minimum_spanning_tree

def sky_figures(unit, mag, uid, sep_deg=20.0, cap=12):
    """unit (N,3) unit vectors, mag (N,), uid (N,) ints -> set of (uid_a, uid_b) edges."""
    order = np.lexsort((uid, mag))                      # brightest first, ties by uid
    u, ids = unit[order], uid[order]
    cos_sep = np.cos(np.radians(sep_deg)); seeds = []
    for i in range(len(u)):
        if all(u[i] @ u[s] <= cos_sep for s in seeds):
            seeds.append(i)
    label = np.argmax(u @ u[seeds].T, axis=1)           # nearest seed by angle
    edges = set()
    for k in range(len(seeds)):
        idx = np.nonzero(label == k)[0][:cap]           # already brightest first
        if len(idx) < 3:
            continue
        sub = u[idx]
        pairs = {(0, 1), (1, 2), (0, 2)} if len(idx) == 3 else {
            tuple(sorted((s[a], s[b]))) for s in ConvexHull(sub).simplices for a, b in ((0, 1), (1, 2), (0, 2))}
        r, c = zip(*pairs)
        w = np.degrees(np.arccos(np.clip((sub[list(r)] * sub[list(c)]).sum(1), -1, 1)))
        tree = minimum_spanning_tree(coo_matrix((w, (r, c)), shape=(len(idx),) * 2).tocsr()).tocoo()
        edges |= {tuple(sorted((int(ids[idx[a]]), int(ids[idx[b]])))) for a, b in zip(tree.row, tree.col)}
    return edges
```

Run over HYG's 600 brightest stars as seen from Earth this gives 56 figures covering 508 stars with median segment 5.8 degrees [C]. Delaunay, MST and RNG for 100 stars take 0.03 s; a few hundred are milliseconds.

### 6.3 Stability and determinism

Moving the observer 20 ly changes 6 of the top 100 stars and keeps 45% of the MST edges (Jaccard on common stars; 29% for the seeded method); at 100 ly 59 stars remain and 24 to 32% of the edges; at 326 ly 38 stars and about 18% [C]. Constellations therefore cannot be galaxy-wide; they belong to a sky, which is a system (a star's planets and moons are within 0.005 pc and share a sky). Store them with the system, as VIEW.3 says, so a later change in the star set does not move them.

Inputs are `(galaxy seed, system uid, sorted visible star uids, algorithm version)`, and every step is a pure function of them, so the result is reproducible under the GEN.39 and GEN.56 rules: no stdlib `random`, and any random choice (such as the optional loop edge) uses `util/draw` seeded by `SHA-256(galaxy seed || "constellation:" system uid || index)` as `galaxy/seed.py` does. Record the algorithm version in the stored rows or the version key. The input list depends on which tiers ran, so build constellations only after the tiers for the chosen `N` have run.

### 6.4 Naming (VIEW.4)

- The codec needs a 19-hex-digit ID. `names/object_id.py::KIND_CODES` stops at 11 (0 and 12 to 63 unused), so `codec_kind_of` cannot recognise a constellation from its bits. Add `constellation: 12` (and, for other galaxies, `galaxy: 13`) and build the ID as the type bits plus the first 70 bits of `SHA-256(galaxy seed || "constellation:" system uid || index)`.
- The codec is injective within a domain (exhaustive for words of up to five digits and 100,000 random 19-digit IDs, [object-ids.md](object-ids.md)), so names are unique galaxy-wide by construction. This answers the open VIEW.4 question "unique per planet or galaxy-wide": galaxy-wide, with no registry.
- The original ask (names from the world's languages and sky cultures) is superseded by the 2026-10-07 plan. If it returns, note the licence: Stellarium's `modern` sky culture says "Text and data: CC BY-SA 4.0" [S, `skycultures/modern/description.md`], and each sky culture declares its own licence [S]; share-alike would apply to a derived list. Codec names need none of this.
- **Bayer designations** [S, Wikipedia "Bayer designation" via search]: a Greek then Latin letter plus the genitive of the constellation name; letters roughly follow brightness with many exceptions; Greek runs out after 24, then Latin letters; Flamsteed numbers follow right ascension. For fictional skies the rule can be exact: sort a constellation's stars by apparent magnitude (ties by uid) and assign alpha, beta, ... omega, then a, b, c. There is no Latin genitive, so write "Alpha <name>" in English form. A star keeps its own name in the system pages; the Bayer name is a sky-chart label, derived and not stored.
- Real figures have more loops and fewer 2-star scraps than any of these methods; a template shape (chain, Y, ring) could be a later polish.

## 7. Rendering plan

- **numpy, scipy and a `zlib` PNG writer (about 15 lines): the first step.** No new dependency. A 1600x800 Hammer PNG is 1.1 MB, a 1000 px stereographic 0.66 MB (noisy images compress poorly). No text.
- **Pillow: add if lines and text must be burned in.** Already installed as scikit-image's dependency (`requirements.lock`: 11.3.0 below Python 3.10, 12.3.0 above; MIT-CMU licence [S, PyPI]). No repo code imports it and `setup.py` does not list it, so declare `pillow` there before importing it. Bundle a font.
- **Client canvas or WebGL with a JSON of stars: the second step.** About 9,000 rows plus the band image gives zoom, hover names, constellation toggles and projection switches without a round trip. The project already does map maths in JavaScript (`static/galaxyprisms.js` mirrors the Python).
- **Not recommended:** cairo (`pycairo`, LGPL-2.1 or MPL-1.1; needs the system library on every platform); matplotlib (heavy, test-only extra); d3-celestial (BSD-3-Clause, last release 0.7.35 on 2020-11-05, depends on d3 ^3.5.17, built for Earth's real sky [S]); Stellarium Web Engine (WebGL, README mentions "Gaia stars database access" [S], licence probably AGPL-3.0 [R], contributors sign a CLA).

Star drawing: peak intensity `min(1, 10^(-0.4 (m - 4.5)))`, core Gaussian sigma 0.7 px plus 0.45 px per magnitude brighter than 4.5, so Sirius (m = -1.46) has sigma 3.4 px and a star at the 6.5 limit peaks at 0.2. These constants are chosen by eye (Stellarium exposes a similar "relative" and "absolute" scale pair [S]; its formula was not confirmed). Splat one point per star with `np.bincount`, blur each channel once with `scipy.ndimage.gaussian_filter`, add the band, and gamma-encode after summing. 1,000,000 faint stars at 1600x800 took 3.3 s; 119,000 stars plus the band lookup 2.7 s [C]. Per-request live rendering is not needed if results are cached by (system uid, epoch bucket, projection, size, limiting magnitude, version key).

Recommendation: stage 1 is a numpy PNG for the all-sky chart; constellation lines and labels are an SVG overlay (the PNG stays cacheable and the labels sharp, and the site already draws maps as SVG); stage 2 is a canvas viewer fed by the same star list.

## 8. Staged plan

1. **Photometry library** (`planetgen/physics/photometry.py`): BC table, `apparent_magnitude`, blackbody colour table, extinction law; tests against the numbers above (Sun -0.07, A_B/A_V 1.3226, Sun B-V to 5,778 K). Pure array functions, no database.
2. **Sky catalogue query** (`planetgen/galaxy/sky.py`): for an observer (system id or a position), arrays of position, L, T and uid from `bright_stars`, generated stars and a backfill to the sky tiers, plus the vectors and light-time step. Needs the tier table and the backfill decision (3.2).
3. **Band** (`galaxy/sky_band.py`): the ray-march on the density, cached by quantised position, plus the dust model; compare against 3.1.
4. **Whole-sky PNG** (Hammer-Aitoff) from stages 2 and 3, written with numpy and `zlib`; `GET /api/systems/<id>/sky.png?proj=hammer&mag=6.5`, cached. This is VIEW.2 plus the picture of VIEW.3 without a horizon.
5. **Constellations** (partition algorithm), stored per system, named by the codec; SVG overlay; Bayer labels.
6. **Horizon view and time of day**, after GEN.104: orientation, twilight, atmosphere. Nebulae and neighbour galaxies as extended sprites.
7. **Canvas viewer** (optional).

## 9. Real catalogues (calibration and tests, not the sky)

- HYG v4.1/4.2: 119,626 rows, 34 MB CSV; CC BY-SA 4.0 from v4.0, CC BY-SA 2.5 earlier [S, astronexus.com/hyg]; share-alike, so a test oracle only, not shipped in the repo.
- Yale Bright Star Catalogue: 9,110 objects, 9,095 stars, V <= 6.5 [S]; the count to match.
- Hipparcos 118,218 stars [S, ESA], used through HYG. Gaia DR3 (over a billion sources, free with credit to ESA/Gaia/DPAC [S], exact licence name unconfirmed) is irrelevant to a generated galaxy.
- Stellarium `modern` sky culture, 88 constellations, CC BY-SA 4.0 and Free Art License [S]: real figures for comparison only.

Calibrate the generator against them (does a generated sky have about 9,000 stars at 6.5? what share are red giants?). Tests can read HYG from a developer machine without committing it.

## Evidence notes

[R] items to verify when paper or catalogue access is allowed (the research environment could read only search-result text and raw GitHub files):

- The BC_V table is a blackbody integral; compare with Flower 1996 or Torres 2010 (hot stars are probably 0.5 mag off).
- The Cardelli, Clayton and Mathis (1989) and Wyman et al. (2013) coefficients were recalled and validated numerically (A_B/A_V = 1.3245 against 1.3226; colour-matching error at most 0.024 against peak 1.78), not read from the papers.
- Dust scale length 2.26 kpc and height 134 pc (Drimmel and Spergel 2001); the 1.0 mag/kpc plane default; the 23 mag/arcsec^2 naked-eye limit; k = 0.2 mag per airmass; the real Milky Way band at 21 to 23 mag/arcsec^2; the local V luminosity density 0.05 L_sun/pc^3.
- Stellarium Web Engine's licence (probably AGPL-3.0) and Gaia DR3's exact licence name (a third-party page says CC BY-NC 3.0 IGO) are unconfirmed.
- The ICRS to galactic matrix used for the demo images was recalled and matches astropy to 1.2e-7.
- Backfill cost per cell is measured (see section 3.2); the cost of the middle 100 to 300 L_sun bands for a whole sky is not.
- HYG is incomplete fainter than about V = 9, so "gained" star counts in the 20 ly parallax test are lower bounds; the "lost" counts and the shifts are unaffected. The unit of HYG's `vx/vy/vz` (pc/yr) was assumed; the median 29 km/s it gave is the only check.
- The kinetics-document matrix check (section 2.6) is [C]; the IAU 2009 report (Archinal et al. 2011) itself was not read.

## Sources

URLs actually used (via search results or raw GitHub fetches):

- Hipparcos count: https://www.esa.int/Science_Exploration/Space_Science/Hipparcos_overview
- Gaia licence and credit: https://cosmos.esa.int/web/gaia-users/license ; https://archives.esac.esa.int/doi/html/data/astronomy/gaia/DR3.html
- HYG database: https://astronexus.com/hyg ; `hyg/CURRENT/hygdata_v41.csv` from https://raw.githubusercontent.com/astronexus/HYG-Database/master/
- Yale Bright Star Catalogue count: https://skyandtelescope.org/astronomy-resources/how-many-stars-night-sky-09172014/ ; https://earthsky.org/astronomy-essentials/how-many-stars-could-you-see-on-a-clear-moonless-night/
- Extinction rates: https://en.wikipedia.org/wiki/Extinction_(astronomy) (as reproduced at https://www.encyclopedia.pub/entry/29910)
- Wyman, Sloan, Shirley 2013: https://research.nvidia.com/labs/rtr/publication/wyman2013simple/ ; https://www.haralick.org/DV/XYZJCGT.pdf
- Ballesteros formula implementations: https://pyastronomy.readthedocs.io/en/latest/pyaslDoc/aslDoc/aslExt_1Doc/ramirez2005.html ; https://en.wikipedia.org/wiki/Color_index
- Star colour perception: https://www.atnf.csiro.au/?p=17326 ; https://www.gresham.ac.uk/node/11799 ; https://www.astronomy.com/observing/observe-the-skys-colorful-stars/
- Bayer and Flamsteed designations: https://www.wikipedia.com/wiki/Bayer_letter ; https://spider.seds.org/spider/Misc/naming.html
- Constellation graph building blocks: https://www.redblobgames.com/x/1812-galaxy-generation/ ; https://portneuf.cose.isu.edu/courses/6673/papers/sewell.pdf (not read)
- Stellarium: https://stellarium.org/doc/head/classStelSkyCulture.html ; the `modern` sky culture on https://raw.githubusercontent.com/stellarium/stellarium/master/skycultures/modern/ ; https://stellarium.org/doc/1.x/classStelSkyDrawer.html ; Web Engine README https://raw.githubusercontent.com/Stellarium/stellarium-web-engine/master/README.md
- d3-celestial: https://raw.githubusercontent.com/ofrohn/d3-celestial/master/LICENSE ; https://registry.npmjs.org/d3-celestial ; constellation lines https://raw.githubusercontent.com/ofrohn/d3-celestial/master/data/constellations.lines.json
- Package pages (PyPI JSON): pillow 12.3.0 (MIT-CMU; 11.3.0 is the last for Python 3.9), numpy, scipy, astropy, pycairo 1.29.2, matplotlib, skyfield 1.55 (MIT), healpy 1.20.1 (GPL-2.0-only)
- Repo files read: `docs/TODO.md`, `docs/design/` (galaxy-coordinate-system, navigation-frames, reproducible-galaxies, object-ids, galaxy-disk-density, galactic-potential, orbital-updates, orbital-solvers-and-integrators, Observational Kinetics for Rotational Vectors), `docs/database-schema.md`, `db/schema.sql`, `galaxy/density.py`, `generation/bright_stars.py`, `tuning.py`, `names/naming_key.py`, `names/object_id.py`, `physics/light_travel.py`, `physics/planets.py`, `web/maps/starmap.py`, `setup.py`, `requirements.lock`, `/mnt/project-files/bright-star-timing/report.md`
