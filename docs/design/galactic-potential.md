# The galaxy's own gravity and galactic orbits

The smooth background pull of the galaxy (bulge, disk and halo) that every
star, remnant and phenomenon feels on its galactic orbit: the model, how to
code it, how to scale it to a generated galaxy, how to step orbits in it and
how to seed velocities so a disk stays a disk. It carries the model and
parameters from Boss's "Computational Astrodynamics.md" (section "Galactic
Potential Modeling and Rotation Curve Equilibrium"), the numbers recomputed
against it, and the places where the source documents and the code disagree.

Informs: GEN.115, GEN.108, GEN.106, GEN.125, GEN.6 (done), GEN.9 (context only)

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

See also: [orbital-updates.md](orbital-updates.md) (the update run, thresholds
and the other solvers), [galaxy-coordinate-system.md](galaxy-coordinate-system.md)
(frame, grid, galactic orbits since 7.37.0) and
[galaxy-disk-density.md](galaxy-disk-density.md) (the star-count model, which
is separate from the gravity model).

Evidence tags: [S] seen in a search result, [C] computed or read from an
installed package during this research, [R] recalled and unconfirmed (see
Evidence notes). The research environment could only read search-result text,
not the papers themselves.

## Decisions already taken

- Boss (2026-10-07 12:25Z): "We'll have to add a galactic gravitational
  gradient but we need to make sure that it's consistent with actual
  science." His research ("Computational Astrodynamics.md") gives the model:
  Hernquist bulge, Miyamoto-Nagai disk, NFW halo (a flattened log halo as a
  setting), nearby point masses as perturbations on top, Plummer softening of
  1 pc.
- Boss (2026-10-07 17:11Z) approved scaling the potential with the galaxy's
  shape, and replaced the document's 1 in-game year step with real time
  ("once set up and configured we follow orbital paths in real time ...
  Default is 1 day = 1 day"), with an option for more time in one go.
- Boss (2026-10-03 05:38Z, GEN.106): a star-like object is updated only when
  it has moved at least 0.01 mpc.
- GEN.115's test targets, from his document: 229.3 km/s at 8.128 kpc and
  213 to 231 km/s from 4 to 20 kpc.

## 1. Findings in brief

- The targets hold, with thin margins (0.47 and 0.24 km/s to the band), and
  the numpy force matches gala to 2e-8 relative [C] (section 2).
- The source document's alternative log halo fails the same targets
  (section 3).
- The existing rotation curve in `galactic_orbit.py` must be replaced in the
  same change: it gives 206 km/s at the Sun and 122 km/s at 2 kpc. Stars
  seeded on it and moved in the new potential become eccentric (Sol's star
  swings between 6.4 and 7.9 kpc; a star at 2 kpc plunges to 1 kpc) [C].
- About 25 lines of numpy, no new library (sections 4 and 5).
- Masses must scale with lengths (section 6).
- Kick-drift-kick leapfrog (the same scheme as Velocity Verlet) with a
  per-star step; no epicycle model (section 7).
- Seed with the asymmetric drift or the disk rearranges within 250 Myr
  (section 8).
- Three Sun distances and stale disk constants need fixing (section 11).

## 2. The model

Source: "Computational Astrodynamics.md", "Galactic Potential Modeling and
Rotation Curve Equilibrium". Axisymmetric, R = sqrt(x^2 + y^2),
r = sqrt(R^2 + z^2), G = 4.300917e-6 kpc (km/s)^2 / Msun. The whole potential
works in kpc, km/s, Msun and time in kpc/(km/s) = 977.79 Myr, converting at
the module edge.

| Component | Model | Potential | Mass | Scale |
|---|---|---|---|---|
| Bulge | Hernquist | `-G M_b / (r + c_b)` | M_b = 1.0e10 | c_b = 0.50 kpc |
| Disk | Miyamoto-Nagai | `-G M_d / sqrt(R^2 + (a_d + sqrt(z^2 + b_d^2))^2)` | M_d = 6.8e10 | a_d = 3.50 kpc, b_d = 0.30 kpc |
| Halo | NFW | `-(G M_h / r) ln(1 + r / r_h)` | M_h = 5.4e11 | r_h = 16.0 kpc |
| Halo (alternative) | Flattened log | `0.5 v0^2 ln(R_c^2 + R^2 + z^2/q^2)` | v0 = 175 km/s | R_c = 2.5 kpc, q = 0.9 |
| Black hole | Plummer-softened point mass | `a = G M r / (r^2 + eps^2)^(3/2)` | 4.3e6 | eps = 1 pc |

Midplane circular speed, `v_c = sqrt(v_b^2 + v_d^2 + v_h^2)`:
`v_b^2 = G M_b R / (R + c_b)^2`, `v_d^2 = G M_d R^2 / (R^2 + (a_d + b_d)^2)^(3/2)`,
`v_h^2 = (G M_h / R) [ln(1 + R/r_h) - (R/r_h) / (1 + R/r_h)]`.

M_h is the NFW scale mass in `Phi = -G M ln(1 + r/r_s) / r`, not a virial
mass. The halo's M200 is 9.0e11 Msun at r200 = 199 kpc (concentration 12.4,
H0 = 70), and the mass inside 100 kpc is 6.8e11 Msun, the same as gala's
`MilkyWayPotential` (6.9e11) [C].

### Rotation-curve check

| Check | Result | Target |
|---|---|---|
| v_c at 8.128 kpc | 229.27 km/s (bulge 68.5, disk 163.6, halo 145.3) | 229.3 |
| Min and max over 4 to 20 kpc | 213.47 (at 20) and 230.76 (at 6.39) | 213 to 231 |
| Local density at 8.128 kpc, z = 0 | 0.1008 Msun/pc3 (disk 0.0914, halo 0.0091, bulge 0.0002) | 0.08 to 0.10 [R] |
| Kz(1.1 kpc) / 2 pi G | 74.9 Msun/pc2 | 70 to 74 [R] |
| Periods at 8.128 kpc: azimuthal, radial, vertical | 218, 158, 83 Myr | about 150 to 170 and 80 to 90 [R] |
| kappa/Omega, nu/Omega | 1.379, 2.63 | a flat curve gives 1.414 |
| Escape speed at 8.128 kpc | 557 km/s | 500 to 600 [R] |

All [C]. The curve v_c(R) in km/s:

| R (kpc) | 2 | 4 | 6 | 8.128 | 10 | 12 | 15 | 20 | 30 |
|---|---|---|---|---|---|---|---|---|---|
| v_c | 190.5 | 223.2 | 230.6 | 229.3 | 226.3 | 223.1 | 218.9 | 213.5 | 205.4 |

It peaks near 6 kpc and falls about 1.2 km/s per kpc from 8 to 20 kpc.

### Reading the parameters

1. A parameter change breaks "213 to 231" (margins 0.24 and 0.47 km/s). G to
   a different digit (4.30091727e-6) changes nothing visible.
2. The disk mass 6.8e10 equals MWPotential2014's and gala's older model. It
   exceeds the stellar mass in `galaxy-disk-density.md` (about 5e10 for disk
   plus bulge, BHG16) because a Miyamoto-Nagai disk stands in for gas, the
   thick disk and the exponential profile's different shape. With b = 0.3 it
   matches an exponential disk of scale length 2.6 kpc best at a = 3.0 to 3.2
   kpc and 1.26 to 1.32 times that disk's mass [C]. So a_d = 3.5 is 1.35 times
   the default scale length, a ratio the scaling rule keeps.
3. The gravity model is not derived from star counts. The bulge here is
   1.0e10 Msun; the density model's bulge holds 31% of the stars (1.4 to
   1.7e10 of 5e10 Msun, BHG16). Nothing couples the two and nothing needs to.

### Against real Milky Way models

Computed by running galpy 1.12.0 and gala 1.12.0 model definitions in a scratch
venv [C]. Irrgang and McMillan are what galpy ships, read from its source.

| Model | v_c at its R0 (km/s) | Local density (Msun/pc3) | Kz(1.1 kpc) / 2 pi G (Msun/pc2) | v_c at 2 / 6 / 15 / 20 kpc |
|---|---|---|---|---|
| Project | 229.3 at 8.128 | 0.1008 | 74.9 | 190 / 231 / 219 / 213 |
| MWPotential2014 (Bovy 2015) | 220 at 8.0 | 0.101 | 71.7 | 187 / 225 / 204 / 198 |
| gala `MilkyWayPotential` | 231.5 at 8.122 | 0.099 | 72.4 | 192 / 235 / 220 / 214 |
| gala `MilkyWayPotential2022` | 229.4 at 8.122 | 0.131 | 71.2 | 180 / 229 / 218 / 212 |
| McMillan 2017 | 233.1 at 8.21 | 0.114 | 74.1 | 196 / 230 / 227 / 223 |
| Irrgang et al. 2013 Model I | 242 at 8.4 | 0.102 | 75.4 | 217 / 243 / 231 / 225 |
| Cautun et al. 2020 | 229 at 8.122 | 0.113 | 73.5 | 206 / 230 / 217 / 208 |

The project's curve sits inside the spread of published models. Gaia-based
measurements give about 229 km/s at 8.1 kpc and a slope near -1.7 km/s per
kpc over 5 to 25 kpc (Eilers et al. 2019) [R], so 209 km/s at 20 kpc against
the model's 213.5; a steeper fall beyond 20 kpc [R] does not matter for a
galaxy that ends near 15 kpc.

### Observed solar values

R0 = 8.122 +- 0.031 kpc is astropy's `Galactocentric` default (reference A&A 615 L15,
GRAVITY 2018) [C]; GRAVITY 2019 and 2021 give 8.178 and 8.275 kpc [R]. Astropy's
v4.0 defaults also give the Sun's velocity (12.9, 245.6, 7.78) km/s and z_sun =
20.8 pc [C]; 245.6 / 8.122 = 30.24 km/s/kpc, matching the Sgr A* proper motion
(about 6.4 mas/yr [R]). Solar peculiar motion is (11.1, 12.24, 7.25) km/s
(Schonrich, Binney, Dehnen 2010) [R]. McMillan 2017 fits R0 = 8.21 kpc and
v0 = 233.1 km/s [S]. Oort constants: A = 15.3, B = -11.9 km/s/kpc (Bovy 2017)
[R], against the project's 14.8 and -13.4 [C].

Two consequences. The Sun's angular speed is about 30.3 km/s/kpc, a galactic
year of about 203 Myr; the project's circular period at 8.128 kpc is 218 Myr
[C]. The "220 to 250 Myr" figure belongs to 220 km/s at 8.5 kpc. And the Sun is
not on a circular orbit (233 to 246 km/s, not 229), which matters only if Sol
is placed from real data; the generated galaxy has no real Sun.

## 3. The log halo does not meet the targets

The source document's alternative halo (v0 = 175 km/s, R_c = 2.5 kpc,
q = 0.9) with its bulge and disk gives v_c of 201, 242, 248, 244, 238, 232,
224, 215 and 203 km/s at 2, 4, 6, 8.128, 10, 12, 15, 20 and 30 kpc [C]. That
is 243.8 km/s at the Sun and 203 to 248 km/s over 4 to 20 kpc, outside both
targets. The flattening q = 0.9 changes nothing in the midplane (the halo's
midplane force depends only on R). A least-squares fit to the NFW total curve
over 2 to 30 kpc gives v0 = 177.2 km/s and R_c = 5.12 kpc: 232.2 km/s at
8.128 kpc, band 213.0 to 232.7, rms 2.6 km/s [C]. It still misses 229.3 by 3
km/s, so a log halo needs its own tuned parameters (or a lower disk mass).
The log potential is also untruncated: no finite escape speed, so stars flung
out never leave.

Recommendation: a halo-type setting, NFW the default and the only type the
229.3 test covers. If the log halo is kept, ship it with v0 = 177 km/s and
R_c = 5.1 kpc.

## 4. Forces, pitfalls, softening, speed

### Formulas

With s = sqrt(z^2 + b_d^2) and D = R^2 + (a_d + s)^2:

- Hernquist: `a = -G M (r + c)^-2 r_hat`. Finite at the centre (G M / c^2 =
  1.7e5 (km/s)^2/kpc) but the direction flips there: a cusp, not a singularity.
- Miyamoto-Nagai: `a_R = -G M R D^(-3/2)`, `a_z = -G M z (a_d + s) / (s D^(3/2))`.
  Smooth for b > 0.
- NFW: `a = -(G M / r^3) [ln(1 + u) - u/(1 + u)] r_vec`, `u = r / r_h`. Tends
  to the finite G M / (2 r_h^2) = 4.5e3 (km/s)^2/kpc at the centre, again with
  a direction flip.
- Flattened log: `a = -v0^2 (x, y, z/q^2) / (R_c^2 + R^2 + z^2/q^2)`; checked
  against gala's `LogarithmicPotential` to 4e-8.

### Numerical pitfalls (tested)

| Pitfall | Where | What happens | Guard |
|---|---|---|---|
| `ln(1+u) - u/(1+u)` cancellation | NFW, small u | float64 relative error 1.5e-12 at 1 pc (harmless); float32 2e-3 at 1 pc, 16% at 0.01 pc | `log1p`; below u = 0.01 use the series `u^2 (1/2 - 2u/3 + 3u^2/4 - 4u^3/5 + 5u^4/6)`; float64 throughout |
| Division by r at r = 0 | Hernquist, NFW | 0/0 gives NaN for a body at the origin | `rs = max(r, 1e-12)` and multiply the factor by the coordinate, so the result is 0 there |
| MN vertical term with b = 0 | disk | `z/s` is 0/0 at z = 0; at b = 1e-9 the force at z = 1e-4 kpc is 1.5e3 (km/s)^2/kpc, a Kuzmin knife edge | `s = max(sqrt(z^2 + b^2), 1e-9)`; reject b below about 0.1 times the grid's layer height |
| Single precision | everything | float32 at 8 kpc resolves 1 pc, 100,000 times the 0.01 mpc threshold | float64 for stored positions and force; float32 only in a throwaway density calculation |

### Softening and the centre (GEN.108)

A Plummer-softened black hole with eps = 1 pc and M = 4.3e6 Msun (Sgr A*'s
order; 4.1e6 to 4.3e6 [R]) has a circular speed that peaks near 81 km/s at 1
pc, then 72 km/s at 3 pc and 43 km/s at 10 pc; unsoftened it would be 136 km/s
at 1 pc and diverge below [C]. Its pull equals the bulge's enclosed-mass pull
inside about 10 pc (bulge M(<r) = M_b (r/(r + c))^2 = 4.3e6 Msun at 10.6 pc),
so the softening touches only the innermost part. The dynamical time at the
softening length is 2 pi sqrt(eps^3 / GM) = 45,000 years, so a star passing
within a few parsecs needs steps of tens of years; the step criterion in 7.4
handles that. A body that really reaches the centre is rare: cap the
sub-steps per update and flag the star when the cap is hit.

Conclusion for GEN.108: the galactic-centre guard is the Plummer softening,
the r = 0 and b = 0 guards above and a sub-step cap. No infinity is left.

### Speed

For 1e6 stars per call, single-threaded numpy, float64 [C]: the straightforward
version takes 0.53 s (median 0.80 s), the tuned one below 0.35 s (median 0.40
s), and float32 0.10 s (error 9e-7, not recommended). A kick-drift-kick step
needs one new acceleration (the last force is reused for the next first
half-kick), so a million stars cost about 0.4 s per step plus array updates;
chunk by 250,000 (about 100 MB with temporaries). In real-time operation only
the due stars move (7.2), so a run costs well under a second per million stars
due; the database read and write is the real cost.

### Tested implementation

Matches the reference to 2e-10 away from the origin and is finite at r = 0,
z = 0 and b_d = 0 [C]. `circular_speed` was added and re-checked at 229.27
km/s.

```python
import numpy as np
G = 4.300917e-6                      # kpc (km/s)^2 / Msun

def galactic_acceleration(x, y, z, p):
    """Hernquist + Miyamoto-Nagai + NFW (+ softened point mass), kpc and km/s -> (km/s)^2/kpc.
    p: Mb, cb, Md, ad, bd, Mh, rh, optional Mbh, eps_bh."""
    R2 = x * x + y * y
    r2 = R2 + z * z
    r = np.sqrt(r2)
    rs = np.maximum(r, 1e-12)
    f = -G * p["Mb"] / ((rs + p["cb"]) ** 2 * rs)            # Hernquist, times (x, y, z)
    s = np.maximum(np.sqrt(z * z + p["bd"] ** 2), 1e-9)
    d = R2 + (p["ad"] + s) ** 2
    d = d * np.sqrt(d)
    fd = -G * p["Md"] / d                                     # MN, times (x, y); z gets (ad + s)/s
    u = rs / p["rh"]
    small = u < 0.01
    us = np.where(small, 1.0, u)
    g = np.where(small, u * u * (0.5 - 2*u/3 + 0.75*u*u - 0.8*u**3 + (5/6)*u**4),
                 np.log1p(us) - us / (1.0 + us))
    fh = -G * p["Mh"] * g / rs ** 3                          # NFW, times (x, y, z)
    if p.get("Mbh"):
        f = f - G * p["Mbh"] / (r2 + p["eps_bh"] ** 2) ** 1.5
    return ((f + fd + fh) * x, (f + fd + fh) * y, (f + fd * (p["ad"] + s) / s + fh) * z)

def circular_speed(R, p):            # midplane, km/s; seeding and stepping share this one function
    ax, _, _ = galactic_acceleration(R, 0.0 * R, 0.0 * R, p)
    return np.sqrt(-ax * R)
```

## 5. Libraries against 25 lines of numpy

PyPI facts read 2026-10-09 [C]. galpy 1.12.0 (BSD) needs Python 3.10+ (1.11.1
is the last for 3.9) and has matplotlib as a hard dependency plus a C
extension. gala 1.12.0 needs Python 3.12+, numpy 2.2+ and astropy 7+, ships
wheels for Linux x86_64 and macOS arm64 only and otherwise needs a compiler. agama 1.0.0 is sdist-only, needs a
C++ compiler and builds action-based models, the wrong tool. Astropy 8.0.1 is
already a dependency (locks use 6.0.1 on Python 3.9 and 6.1.7 on 3.10).

The project supports Python 3.9+ on Linux and macOS (`setup.py`).
galpy would add matplotlib and a compiled extension for three analytic
potentials; gala would break Python 3.9 to 3.11. Neither gives
anything the numpy version lacks (2e-8 against gala, 0.4 s per million stars).
Recommendation: no runtime dependency; an optional test against gala when it
is installed, skipped otherwise.

Astropy's `Galactocentric` frame (for VIEW.2, not for the potential) puts the
Sun at negative x; a 200,000-point conversion to ICRS took 0.21 s [C]. With the
Sun at x = -R0 moving along +y, the real Milky Way's angular momentum along z is
negative: it rotates clockwise seen from the north galactic pole [C,
arithmetic]. The project's galaxy rotates counterclockwise seen from +Z
(`system_position.galactic_velocity_ms`, `galaxy-coordinate-system.md`
section 1), so its "galactic north" is the real galaxy's south. Mapping to
real sky coordinates needs a mirror in z or y.

## 6. Scaling the potential to the generated galaxy

Code defaults (`web/generate_page.py`, `galaxy/density.py`): disk scale length
2,600 pc, scale height 300 pc, bulge scale radius 1,580 pc (the Dwek bar
scale) and an edge near 50,000 ly (15.3 kpc; the outline reaches 50,300 ly).
The potential's lengths were chosen with the Milky Way in mind.

With s = disk_scale_length / 2,600 pc, scaling every length by s and keeping
the masses (the rule now written in `orbital-updates.md` and GEN.115) keeps the
curve's shape but not its speed: v_c at the same R/s goes as s^-1/2 [C].

| s | Disk scale length | v_c at 3.15 L, masses fixed (km/s) | Masses x s | Masses x s^2 |
|---|---|---|---|---|
| 0.02 | 52 pc | 1,621 | 229.3 | 32.4 |
| 0.1 | 260 pc | 725 | 229.3 | 72.5 |
| 1 | 2.6 kpc | 229.3 | 229.3 | 229.3 |
| 4 | 10.4 kpc | 115 | 229.3 | 458 |
| 10 | 26 kpc | 72.5 | 229.3 | 725 |

Fixed masses give 1,600 km/s for a tiny galaxy and a sluggish 72 km/s for a
large one. The curve's shape holds in every column (the 4 to 20 kpc band,
rescaled to 4s to 20s kpc, spans 0.93 to 1.01 times the value at the Sun).

### Recommended rule

1. **Lengths.** a_d = 3.5 kpc x s_d and r_h = 16 kpc x s_d, with
   s_d = disk_scale_length / 2,600 pc. b_d comes from the shape's
   `disk_scale_height_pc` (default 300 pc = 0.3 kpc exactly) and
   c_b = 0.5 kpc x bulge_scale_radius / 1,580 pc. The density and gravity
   models then describe the same disk and bulge. "Everything by disk scale
   length" is the special case at default ratios.
2. **Masses.** M_i = M_i(default) x s_d^p, p = 1 by default (a
   `galactic_mass_exponent` setting, 1 to 2). p = 1 holds 229.3 km/s at 3.15 L
   and the same shape in R/s [C]. p = 2 (constant surface density, Freeman's
   law) gives v proportional to s^0.5 and the Tully-Fisher scaling
   M proportional to v^4 [R, reasoning]; it is more physical but gives 725
   km/s at s = 10 and needs a clamp. A friendlier setting exposes "rotation
   speed at the Sun's position, km/s" (default 229.3), computes the mass scale
   from it and clamps to 30 to 600 km/s. Nothing couples star counts to the
   potential, so only plausible times and speeds are at stake.
3. **Black hole.** Keep 4.3e6 Msun and eps = 1 pc unscaled, but clamp the mass
   below 1% of the bulge mass.
4. **Derived.** The orbital period at the Sun's place scales with s_d for p = 1:
   218 Myr at s = 1, 22 Myr at 0.1, 2.2 Gyr at 10 [C]. The 0.01 mpc threshold is
   in absolute parsecs, so a small galaxy has more due stars per run and a
   huge one fewer; that is correct. At 229 km/s the interval is 15.6 days at
   any size.
5. **Guards.** Reject s_d outside about 0.02 to 20 and b_d / a_d above 0.5
   (the disk becomes a Plummer sphere); log the derived v_c at the calibration
   radius; keep every component speed below 0.01 c (3,000 km/s).

At s_d = 1 and default ratios this reproduces Boss's numbers exactly, so the
GEN.115 tests at the default shape do not change.

The parameters are a function of the shape, so for GEN.9 they belong with the
`galaxy_shape` row, derived at read time rather than stored; a second galaxy
then has its own potential with no schema change. Other galaxies seen from
this one (GEN.9 stage 1) need no potential: they are distant objects with a
position and a catalogued velocity.

## 7. Orbit integration

### 7.1 Frequencies

At 8.128 kpc [C]: Omega = 28.21 km/s/kpc, kappa = 38.89, nu = 74.29 (periods
218, 158 and 83 Myr). At 4 kpc the periods are 110, 73 and 40 Myr; at 15 kpc
421, 311 and 182 Myr. In the bulge they fall to a few Myr (19 Myr at 0.5 kpc).

### 7.2 What a real-time update looks like

A star at 230 km/s moves 0.64 millipc per day, so the 0.01 mpc threshold
(GEN.106) is reached after about 15.6 days; `next_update_due` is that long
after the last update and a daily run touches about 1/15 of the stars. A slow
halo star near apocentre (50 km/s) is due every 72 days.
`next_update_due = now + 0.01 mpc / |v|` is adequate, since speeds change over
Myr; recompute whenever the star moves.

A single leapfrog step of 15.6 days has a relative frequency error
(Omega dt)^2 / 24 of order 1e-11. A catch-up run after a long gap, or the
`--extra-time` option, takes sub-steps by the criterion in 7.4. Do not loop one
day at a time (about 200 microseconds per Python call per step). Float64
roundoff over 1e6 daily steps (2,738 years) for 8 test stars, against an 80-bit
long-double run of the same steps, was 9e-12 pc at both 1e5 and 1e6 steps [C],
so there is no rounding drift to engineer around. Sector-relative milliparsec
storage is fine: an absolute position rebuilt in float64 at 8 kpc resolves
1.8e-12 pc, far below the 1e-5 pc threshold, provided the force is evaluated on
the reconstructed galactic position and the displacement added back.

### 7.3 Leapfrog accuracy

Kick-drift-kick leapfrog in the full model with the softened black hole, 60
thin-disk stars (sigma_R 35, sigma_z 20 km/s, lag 15 km/s), 1 Gyr (about 4.6
turns), against a dt = 0.01 Myr reference that agrees with SciPy DOP853 at
rtol 1e-12 to 0.003 pc [C]:

| dt (Myr) | max dE/E | max dLz/Lz | Position error at 1 Gyr, median / max |
|---|---|---|---|
| 0.5 | 2.2e-5 | 8e-15 | 4 / 15 pc |
| 1 | 8.9e-5 | 7e-15 | 16 / 59 pc |
| 2 | 4.2e-4 | 6e-15 | 64 / 235 pc |
| 5 | 1.7e-3 | 3e-15 | 397 / 1,522 pc |
| 10 | 1.3e-1 | 2e-15 | 1,600 / 8,044 pc |

Lz is conserved to roundoff at any step (the axisymmetric force is central in
the plane). The energy error is bounded, goes as dt^2 and runs away only as dt
nears the shortest period along the orbit. Other populations need much
smaller fixed steps: thick disk at 2 Myr gives dE/E 7.5e-2, halo-like at 0.5
Myr 4.8e-2, inner at 1 Myr 0.62 [C]. A circular star at 8.128 kpc stays
within 0.04 pc of its radius over 5 Gyr at dt = 0.2 Myr, 0.9 pc at 1 Myr, 22 pc
at 5 Myr and 89 pc at 10 Myr [C].

### 7.4 Step criterion

Per-star step, chosen once per update from the star's state:
`dt_i = min(eta * min(t_peri, t_local), dt_max)`, with `t_peri = r_p / v_p` the
pericentre crossing time (r_p from L = r_p v_p and the energy, by a few
fixed-point iterations of `v_p^2 = 2 (E - Phi(r_p))`) and
`t_local = sqrt(r / |a|)`. Round each step down to dt_max / 2^k and integrate
the stars of each k together (block steps). Max dE/E over 1 Gyr [C]:

| eta | Thin disk: median dt, steps | Thin disk dE/E | Thick disk dE/E | Halo-like dE/E | Inner (0.2 to 1.5 kpc) dE/E |
|---|---|---|---|---|---|
| 0.20 | 4.9 Myr, 317 | 8e-4 | 6e-4 | 5.9e-2 | 2.5e-2 |
| 0.10 | 2.4 Myr, 627 | 2e-4 | 1.6e-4 | 1.7e-3 | 4.4e-4 |
| 0.05 | 1.2 Myr, 1,253 | 4.9e-5 | 3.9e-5 | 4.2e-5 | 2.6e-5 |
| 0.02 | 0.5 Myr, 3,253 | 8.6e-6 | 5.2e-6 | 3.6e-6 | 3.1e-6 |

Recommended: eta = 0.05 and dt_max = 5 Myr (about 1,250 steps per Gyr for a
thin-disk star, 12,000 to 23,000 for halo and bulge stars). Position error at
eta = 0.05 is a median 14 pc (thin), 5 pc (thick), 1 pc (halo) and 2.4 pc
(inner), with the worst inner star 637 pc off after a Gyr (near-chaotic
orbits). That is fine for a game, since positions stay on an orbit of the
right energy and angular momentum, but long catch-up runs are not
reproducible against a different eta: fix eta and dt_max as constants and
record them in the run log. Changing the step between updates (not within one)
keeps the scheme symplectic per update; energy wanders by bounded errors of
1e-5 per update at worst, well under 1e-3 after 1,000 updates.

Cost: a one-time Gyr catch-up for a million stars is about 1e9 star-steps at
0.4 microseconds, about 400 s. Cap the sub-steps at about 50,000 per star per
update; a star that would need more (orbiting the black hole at a few parsecs
through a long update) keeps its current state, is flagged and logged. That
is the GEN.108 "limit" for the centre.

### 7.5 Epicycles are not the update

Textbook epicycle formulas (guiding radius from Lz, radial oscillation at
kappa, the azimuthal drift term, vertical oscillation at nu) against the exact
orbit for 80 thin-disk stars [C]: median error 213 pc at 25 Myr, 1.10 kpc at
100 Myr, 3.5 kpc at 500 Myr and 6.7 kpc at 1 Gyr (maxima 2.98, 7.5, 13.1 and
16.8 kpc). Even a star with vR = 5 km/s is 35 pc off after 25 Myr, because the
real frequencies depend on amplitude, and a 300 pc vertical oscillation is
off by 250 pc within 25 Myr because the Miyamoto-Nagai vertical force is far
from linear there. Use epicycles for initial conditions and diagnostics, not
for the update.

One exact shortcut remains: a body on a planar circular orbit (zero peculiar
velocity, z = 0) can be rotated by `Omega(R) t`, which is what the GEN.6 code
does for every star today. It is right for bodies with no random motion and
must give way to integration for stars with dispersion.

## 8. Seeding velocities

### 8.1 Asymmetric drift

A population with radial dispersion sigma_R does not rotate at the circular
speed. From the Jeans equation for an exponential disk (density scale length
L, sigma_R^2 falling as exp(-R/L_s)) with sigma_phi^2 / sigma_R^2 =
kappa^2 / (4 Omega^2):

`v_c^2 - vbar_phi^2 = sigma_R^2 [R/L + R/L_s - 1 + kappa^2/(4 Omega^2)]`,
`v_a = v_c - vbar_phi`.

With L = L_s = 2.6 kpc the bracket at the Sun is 5.73, so v_a is 5, 16, 27 and
50 km/s for sigma_R of 20, 35, 45 and 60 km/s [C]. Stromberg's
`v_a = sigma_R^2 / 80 km/s` gives 5, 15, 25 and 45 [R]. Observed lags are a few
km/s for young stars, 15 to 25 for old thin-disk stars and 40 to 50 for the
thick disk [R]. Use the Jeans form: its bracket grows outward and shrinks
inward.

Without the drift [C]: 60,000 thin-disk test stars from an exponential disk
(L = 2.6 kpc, sigma_R = 35 km/s at the Sun falling as exp(-(R - R0)/2L),
capped at 90, sech^2 vertical profile of 300 pc scale height), evolved 1 Gyr in
the smooth potential. Change in star count by radius bin (R = 2-4, 4-6, 6-8,
8-10, 10-12, 12-14, 14-16 kpc):

| Seed | At 250 Myr | At 1 Gyr |
|---|---|---|
| Mean v_phi = v_c - v_a (Jeans) | -3.9, -3.2, +2.8, -2.4, +7.1, +10.8, +1.8 % | -4.3, -1.6, +0.6, +1.5, +5.0, +7.9, +7.7 % |
| Mean v_phi = v_c (no drift) | -15.1, -4.9, +5.4, +22.7, +77.7, +24.7, +3.9 % | -15.7, -5.2, +12.0, +25.9, +39.8, +29.7, +44.9 % |

Without the drift, mean v_phi at 10 to 14 kpc falls from about 222 to 200
km/s, the radial dispersion there rises from 21 to 37 km/s and the outer disk
fills up. With it the profile changes by under about 10%. Either way the
vertical rms thickness grows by 10 to 20% because the local velocity
distribution is not a function of the actions; that is the limit of a
position-space seed (a quasi-isothermal distribution function in actions,
Binney 2010, galpy's `quasiisothermaldf`, would remove it [R]). Accept the
transient or generate the scale height 15% lower. Do not pre-evolve: it costs
the whole catch-up and is deterministic only with eta frozen.

### 8.2 Dispersions per population

The populations already exist in `galaxy-disk-density.md` section 4 (young,
intermediate, old with exponential scale heights 50, 150 and 300 pc; thick
disk 900 pc; bulge). The vertical dispersion that holds a layer of a given
scale height follows from the potential: matching the second moment of an
isothermal layer `rho(z) ~ exp(-(Phi(z) - Phi(0)) / sigma_z^2)` to a
sech^2(z/2h) layer gives [C]:

| h (pc) | sigma_z at 4 kpc (km/s) | at 8.128 kpc | at 12 kpc |
|---|---|---|---|
| 50 (young) | 13.1 | 6.4 | 3.9 |
| 150 (intermediate) | 31.7 | 15.9 | 9.9 |
| 300 (old) | 49.5 | 26.0 | 16.5 |
| 900 (thick) | 86.5 | 53.4 | 36.5 |

The values at the Sun match the observed age-velocity relation roughly (young
6 to 10, old thin 20 to 25, thick 35 to 45 km/s [R]; the thick value runs 20%
high because the thick disk here has a 900 pc scale height). sigma_z must fall
outward (about as exp(-R / 5 kpc)) to keep a constant scale height; constant
sigma_z flares the disk. Tabulate sigma_z(R, h) from the potential at galaxy
setup (40 by 8 points, once per galaxy).

Radial dispersion: sigma_R = sigma_z / 0.55 for thin-disk populations (0.65 for
old), falling as exp(-(R - R0) / 2L) outward and capped at 90 to 100 km/s
inside; sigma_phi = sigma_R kappa / (2 Omega) = 0.69 sigma_R at the Sun. At the
Sun (sigma_R / sigma_z / lag, km/s) [R for the observed ranges]: young 15 / 6 /
2, intermediate 28 / 16 / 9, old 40 / 26 / 20, thick 60 / 45 / 45; bulge an
isotropic rotator with sigma about 110 and mean rotation from the Jeans
equation; halo without rotation, sigma about 120. Draw velocities as Gaussians;
reject draws that put the apocentre beyond the galaxy edge or the pericentre
inside 100 pc.

## 9. Test values

For the GEN.115 tests at the default shape:

- v_c(8.128 kpc) = 229.3, `approx(abs=0.1)`.
- The 213 to 231 band checked on a 2 kpc grid (4, 6, ..., 20), not a fine
  scan (the continuous maximum is 230.76 and minimum 213.47).
- `abs(dLz/Lz) < 1e-10` at any step size.
- Energy drift below 1e-3 over 1 Gyr at eta = 0.05.
- A circular star at 8.128 kpc keeps its radius to 1 pc over 5 Gyr at dt <= 1
  Myr.
- Kz(1.1 kpc) about 75 Msun/pc2 and local density about 0.10 Msun/pc3, 5%
  tolerance.
- The acceleration at (0, 0, 0) is finite, also with b_d = 0.
- Optional: agreement with gala to 1e-6 relative when gala is installed.
- A thin-disk sample seeded with the drift changes its radial count profile by
  less than about 10% in 250 Myr.

## 10. How this fits planetGen

- `galaxy/galactic_orbit.py::calculate_galactic_orbit(distance_ly)` and the
  constants `GALACTIC_ROTATION_FLAT_VELOCITY_KMS = 220` and
  `GALACTIC_ROTATION_CORE_RADIUS_PC = 3000` in `physics/constants.py` are what
  GEN.115 replaces. Their comment claims "~206 km/s and a ~236 million-year
  orbit" at Sol's 25,800 ly (7.91 kpc); the new potential gives 229.6 km/s and
  212 Myr there [C]. The old curve gives 70 km/s at 1 kpc where the potential
  gives 168. The function is called from `generate_galactic_orbit_fields`,
  i.e. by every star, remnant and phenomenon at generation, so the potential's
  circular speed must feed the same function.
- `db/store.py::advance_galactic_positions` rotates each system, phenomenon
  and stand-alone facility about Z by `_galactic_turn(elapsed_years,
  period_gy)` and rotates the stored velocity by the same angle: the GEN.6
  scheme, exact for a planar circular orbit, which GEN.115 replaces for stars
  with peculiar velocity. Sector refiling (`_SectorIndex.sector_at`) stays;
  only the motion rule changes.
- `star_systems.velocity_*_kms` (schema v61), the `planets`, `moons` and
  `comets` velocity columns (v60) and `facilities.velocity_*_kms` (v64,
  GEN.125) already hold galactic velocities. With the potential a star's
  velocity becomes `v_c(R)` along the tangent plus a peculiar part from the
  dispersions above (runaways keep `runaway_speed_kms`). Existing databases
  hold velocities on the old curve, so a migration or `reseed-velocities` step
  must recompute them or every old star is eccentric. GEN.125's stand-alone
  facilities take `v_c(R)` from the same function with no dispersion; the
  rotation scheme stays exact for them unless they sit above the plane, where
  they oscillate vertically.
- `physics/sector_path.py` (GEN.123) integrates a test particle against the
  sector's point masses with 1 pc softening; adding the galactic acceleration
  keeps a path consistent when it leaves a sector at low speed.
- `galaxy/keepout.py::galactic_hill_radius_km` could use the potential's
  enclosed mass (6.8e11 Msun inside 100 kpc, 2.1e11 inside 20 kpc) instead of
  the single `MILKY_WAY_MASS` (1.15e12 Msun). A possible follow-up.

## 11. Corrections to the source documents

"Computational Astrodynamics.md" is Boss's upload and is not edited. Where it,
`orbital-updates.md`, the TODO text or the code disagree with the checks above,
this document is the reference.

1. **The log halo fails the targets.** The document's table lists the flattened
   log halo (v0 = 175 km/s, R_c = 2.5 kpc, q = 0.9) as an alternative;
   `orbital-updates.md` section 10.1 and GEN.115 repeat it as a setting under the
   same tests. It gives 243.8 km/s at 8.128 kpc and 203 to 248 km/s over 4 to
   20 kpc (section 3) and needs v0 about 177 and R_c about 5.1, or its own tests.
   The document's "flat within a narrow velocity band (213-231 km/s)" holds
   for NFW only.
2. **The scaling rule keeps the shape but not the speed.** `orbital-updates.md`
   section 11 and GEN.115 scale every length by the disk scale length ratio
   and keep the masses. Speed then goes as s^-1/2: 725 km/s at one tenth the
   size, 72 km/s at ten times. Masses must scale with the lengths (default
   exponent 1), and b_d and c_b should follow the shape's own scale height and
   bulge radius (section 6).
3. **Stale constants.** `orbital-updates.md` section 11 quotes disk scale
   length 2,800 pc, bulge scale radius 200 pc and radius 15 kpc, and divides by
   2,800 pc. The code defaults are 2,600 pc, 1,580 pc and about 15.3 kpc
   (`web/generate_page.py`, `galaxy-disk-density.md` section 2); the divisor is
   2,600 pc. The 200 pc bulge was revision 2's, replaced by the COBE bar
   (`galaxy-disk-density.md` section 5).
4. **Three Sun distances.** `physics/constants.py` has
   `GALACTIC_CENTER_DISTANCE_LY = 25800` (7.91 kpc). The source document, GEN.115
   and `orbital-updates.md` use 8.128 kpc (the document allows 8.128 to 8.200).
   `tuning.GALAXY_SOLAR_RADIUS_TO_SCALE_LENGTH = 8200 / 2600` puts the density
   model's Sun at exactly 8.2 kpc, which the docs round to "3.15 L" (that
   rounds to 8.19). No test uses 8.128 yet, because GEN.115 is unbuilt. The
   potential itself does not care; the test calibration radius and Sol's
   placement do. Pick one reference and move the others to it.
5. **"Singularity-free" and "vanishes at the origin".** The document calls the
   Hernquist derivatives singularity-free and says Plummer softening makes the
   net acceleration vanish at the origin. The Hernquist and NFW forces are
   finite at r = 0 but flip direction (a cusp), so the r = 0 and b = 0 guards
   in section 4 are still needed; the zero at the origin belongs to the
   softened black hole alone.
6. **Macro step and integrator.** The document puts host stars on Velocity
   Verlet with a 1 in-game year macro step. Boss (2026-10-07 17:11Z) replaced
   the year with real time. Velocity Verlet and kick-drift-kick are the same
   scheme at a fixed step; section 7.4 adds a per-star step because a fixed step
   fails for eccentric, halo and bulge stars.
7. **The NFW mass is a scale mass.** The "characteristic halo mass parameter"
   5.4e11 Msun corresponds to a virial mass of 9.0e11.
8. **The old rotation-curve comment.** `physics/constants.py` calls ~206 km/s
   and ~236 Myr "close to the real Sun's measured ~220-240 km/s and
   ~225-250 million-year galactic year". The Sun's angular speed from the
   Sgr A* proper motion gives about 203 Myr, and the new model gives 218 Myr
   at 8.128 kpc (section 2).
9. **"Orbits keep their radii".** The document says the potential keeps stellar
   radii over long integrations. That holds only if stars are seeded from the
   same potential with the asymmetric drift: seeding at the circular speed
   shifts outer-disk counts by +25 to +78% within 250 Myr (section 8), and the
   old curve makes every star eccentric (section 1).
10. **Observational claims.** The document attributes 229 to 232 km/s to Gaia
    DR3, LAMOST and DESI; the research could not open those papers (the Gaia
    value, about 229 at 8.1 kpc, is [R]). The tests therefore state "at the
    default shape" and do not claim an observation.
11. **Rotation sense.** `galaxy-coordinate-system.md` section 1 says the frame
    "does not match any real sky frame"; the rotation sense is one more
    difference (counterclockwise from +Z against clockwise from the real north
    pole, section 5).

## Evidence notes

Web search covered three queries before the shared budget ran out, and the
environment could read only search-result text. To verify when paper access is
allowed:

- GRAVITY 2019 R0 = 8.178 kpc and 2021 R0 = 8.275 kpc with uncertainties; Sgr
  A* mass 4.1e6 and 4.3e6 Msun. The 2018 value is confirmed only through
  astropy's reference link.
- Reid and Brunthaler proper motion (6.379 mas/yr in 2004, 6.411 in 2020) and
  the derived 30.2 to 30.4 km/s/kpc (astropy's own defaults give 30.24).
- Schonrich, Binney and Dehnen 2010 solar motion (11.1, 12.24, 7.25).
- Eilers et al. 2019 (229.0 km/s, slope -1.7 km/s/kpc) and the steeper Gaia DR3
  decline beyond 20 kpc.
- Kz(1.1 kpc) of 71 and 74 Msun/pc2 (Kuijken and Gilmore 1991; Holmberg and
  Flynn 2004), local density 0.08 to 0.10 Msun/pc3 (Boss's range), Oort
  constants (Bovy 2017), escape speed 500 to 600 km/s.
- Allen and Santillan 1991 parameters (coded from memory; they reproduce the
  design 220 km/s at 8.5 kpc, which is the only check); left out of the
  comparison table for that reason.
- Observed age-velocity relations, thick-disk kinematics and the Stromberg
  relation (section 8).
- Irrgang 2013 and McMillan 2017 values are what galpy ships, not a re-reading
  of the papers.
- The Tully-Fisher and Freeman's-law reasoning in section 6 is from memory.
- That astropy's latest defaults equal the v4.0 set was read from astropy 8.0.1;
  the locked 6.0.1 on Python 3.9 was not run.

## Sources

- McMillan 2017 (MNRAS 465, 76), search result [S]:
  https://academic.oup.com/mnras/article/465/1/76/2417479 and
  https://ar5iv.labs.arxiv.org/html/1608.00971 (R0 8.21, v0 233.1; abstract
  8.20 +- 0.09 and 232.8 +- 3.0; r_h 19.6 kpc).
- gala docs, search result (structure only; numbers from running the library) [S]:
  https://gala-astro.readthedocs.io/en/latest/api/gala.potential.potential.MilkyWayPotential2022.html and
  https://gala-astro.readthedocs.io/en/latest/api/gala.potential.potential.BovyMWPotential2014.html
- Bovy 2015 galpy paper, search result (parameters from galpy 1.12.0
  `mwpotentials.py`) [S]: https://arxiv.org/pdf/1412.3451
- PyPI JSON [C]: https://pypi.org/pypi/galpy/json,
  https://pypi.org/pypi/gala/json, https://pypi.org/pypi/astropy/json,
  https://pypi.org/pypi/agama/json (and per-version JSON for galpy 1.11.1 and
  1.12.0, gala 1.8.1, 1.10.1 and 1.12.0).
- Installed packages read directly [C]: galpy 1.12.0, gala 1.12.0, astropy 8.0.1
  (`galactocentric_frame_defaults`).
- Project files read: `docs/TODO.md` (GEN.106, 108, 115, 109, 9, VIEW.2),
  `docs/design/orbital-updates.md`, `docs/design/Computational Astrodynamics.md`,
  `docs/design/galaxy-coordinate-system.md`, `docs/design/galaxy-disk-density.md`,
  `src/planetgen/physics/constants.py`, `galaxy/galactic_orbit.py`,
  `galaxy/density.py`, `tuning.py`, `galaxy/system_position.py`, `db/store.py`,
  `web/generate_page.py`, `setup.py`, `requirements.lock`.
- The computation scripts live in the research scratchpad, not in the repo. The
  rotation curve, density, Kz and the code block were rerun while writing this
  document and reproduced (229.27 km/s, band 213.47 to 230.76, 0.1007 Msun/pc3,
  Kz about 74.8 to 74.9, finite acceleration at the origin).
