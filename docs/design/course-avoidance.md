# Course avoidance: keep-out radii and the path-bending algorithm

What radius a course keeps from each kind of object, where the built rule (`galaxy/keepout.py`, NAV.24) has to change, the algorithm that bends a straight course around keep-out spheres (NAV.26), how moving bodies inside a system are handled (NAV.27) and the `adjusted` course object that the NAV page, the maps and saved courses share (NAV.28). It also records which parts of Boss's uploaded astrodynamics documents carry over and which do not. The route itself (which systems are the stops) is in [course-routing.md](course-routing.md); the frames a course is written in are in [navigation-frames.md](navigation-frames.md).

Informs: NAV.6, NAV.22, NAV.24, NAV.25, NAV.26, NAV.27, NAV.28, NAV.51, NAV.5, NAV.21, NAV.4, NAV.17, MAP.121
Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Prototype and test scripts are in the research scratchpad (`nav05/`: planner `avoid.py`, numbers `physics.py`, `collide.py`, `moving.py`, `moving2.py`, `nebfrac.py`, tests `t1.py` to `t20.py`). The scripts are research code outside the repo. The research environment could only read search-result text, not the papers themselves; [R] items are listed under Evidence notes.

## Decisions already taken

- NAV.6 (Boss, 2026-10-01 20:55Z): "we need to factor gravitational bodies into the course. A ship piloting would adjust the course to avoid falling into the gravitational field of objects it knows about. ... attempting to stay out of the Hill sphere of each object. We're also going to use this within the sector and within the star system."
- NAV.51 (Boss, 2026-10-09 01:02Z): "Courses should avoid asteroid fields." Nebulae and supernova remnants were not ruled on; the course passes through them with a note until Boss decides.
- The NAV chain is parked (Boss, 2026-10-09).
- Related, from the route work (Boss, 2026-10-02 01:53Z): no hop limit, and a hop through unknown space is flagged. Avoidance therefore bends each hop of a route one at a time and does not change which systems are the stops.
- Not a Boss decision, a pre-planning default recorded in the NAV.6 text: a lone star or object uses its galactic Hill radius.

## Summary

- **Keep-out rule.** NAV.24 is built and its rule is sound for planets, moons, stars, systems and the compact objects (one default per kind, section 2). The galactic Hill spheres of stars fill only about 6% of the disk volume near the Sun, so they are separate obstacles and not a wall.
- **Four fixes to the built rule** (section 3): asteroid fields need a hard keep-out at `radius_ly` (built code returns none, contradicting NAV.51); the supermassive hole and quasar get the event horizon (0.08 AU) and need the sphere of influence (1.7 pc); a moon's stored `hill_radius_km` is computed against the star's mass (confirmed in the code); the nebula and remnant note gives the wrong reason for pass-through.
- **Algorithm** (section 5): straight line, then the exact one-sphere tangent-and-arc solution wherever a leg crosses a sphere, then polish by block-coordinate descent. In 150+ scenes it was always valid, never longer than a lattice A*, and 10 to 100 times faster (1 to 8 ms for up to 600 spheres).
- **Moving bodies** (section 6): at warp 1 or faster, departure-time positions are enough. **Exposure** (section 7): an `adjusted` object beside `direct` and `route`.

## 1. What is in the repo (checked)

| Piece | What it does | Note |
|---|---|---|
| `galaxy/keepout.py`, `queryDb.keep_out_radius` | One rule per kind, returns `KeepOut(radius_km, basis, note)`; `GET /api/objects/<ref>` shows it as `keep_out` | Asteroid field returns none (needs the NAV.51 change); SMBH and quasar floor is the horizon radius |
| `stars.system_perimeter_km` | `r_gal * (M / (3 M_MW))^(1/3)`, `MILKY_WAY_MASS` = 1.15e12 solar masses (`Star.calculate_system_perimeter`) | Sun at 8 kpc: 0.53 pc = 1.7 ly = 1.09e5 AU [C] |
| `planets.hill_radius_km`, `moons.hill_radius_km` | `a (m / 3M)^(1/3)` (`orbits.calculate_hill_sphere`) | Planet: central mass is the star, correct. Moon: also the star's mass (section 3.3) |
| `asteroid_fields.radius_ly` | Bounding sphere, 0.001 to 1 ly (63 AU to 63,000 AU), `tuning.ASTEROID_FIELD_RADIUS_RANGE_LY` | No number density stored, only `density` = sparse, typical or dense |
| `nebulae.radius_ly`, `nebula_shape.NebulaShape.contains` | Bounding radius 0.1 to 200 ly; irregular metaball shape | A containment test exists, so a chord length can be measured |
| `supernova_remnants.radius_ly` | `0.35 ly * age^0.4` (`SEDOV_TAYLOR_RADIUS_COEFFICIENT_LY`) | Coefficient reviewed in [nebula-and-asteroid-field-classes.md](nebula-and-asteroid-field-classes.md) |
| `physics/body_positions.positions_at(scene, years)` | Every body of a system at any time | Built (MAP.70); NAV.27 uses it as is |
| `galaxy/geometry.sectors_along_segment` (NAV.38) | Sectors a straight segment crosses | Starting point for NAV.25's broad phase |
| `web/nav_page.galaxy_course` | Course for the Galaxy Map: `points` `[{name, url, role, x, y, z}]` in parsecs, `direct`, `readout` | The extension point for the adjusted form |

Sector size is `DEFAULT_SECTOR_EDGE_PC = 4` (13.05 ly), in rings, layers and slots, not true cubes.

## 2. Keep-out radius per kind

One default per kind. Where a kind has two rows (stars, neutron stars) the first applies at galaxy and sector scale and the second inside a system, which is what "within the sector and within the star system" (Boss, NAV.6) needs.

| Kind | Keep-out radius | Class | Data in the DB |
|---|---|---|---|
| Star, white dwarf, O/B giant (galaxy and sector scale) | galactic Hill radius `r_gal (M / 3 M_MW)^(1/3)`: 0.25 to 2.1 pc | Hard | `stars.system_perimeter_km` |
| Binary or multiple system (galaxy and sector scale) | the larger of the pair's combined-mass galactic Hill radius and each star's | Hard | `star_systems.binary_system_perimeter_km`, `stars.system_perimeter_km` (as built) |
| Star inside its own system | `max(10 R_star, sqrt(L / (4 pi F_lim)))`, F_lim = 50 kW/m^2 (tunable) | Hard; exempt when an end lies inside it | `stars.radius_km`, `stars.luminosity_w` |
| Planet | `hill_radius_km` | Hard (in system and sector) | `planets.hill_radius_km` |
| Moon | not needed from outside (the planet's sphere covers it); inside the planet's sphere, Hill radius against the planet's mass | Hard, local only | `moons.distance_km`, `moons.mass_kg`, planet mass (stored `moons.hill_radius_km` is wrong) |
| Asteroid belt | none; note when crossed | Pass-through | `asteroid_belts` limits |
| Asteroid field | `radius_ly` (bounding sphere) | Hard (Boss) | `asteroid_fields.radius_ly` |
| Comet, interstellar comet | none | Pass-through | |
| Rogue planet | galactic Hill radius, floor its own radius | Hard (tiny) | `rogue_planets.mass_kg`, `radius_km`, `galactic_radius_pc` |
| Stellar or intermediate black hole | galactic Hill radius, floor the horizon | Hard | `black_holes.mass_solar`, `galactic_radius_pc` |
| Supermassive black hole | `G M / sigma^2` with the generator's sigma = 100 km/s (1.7 pc at 4e6 solar masses) | Hard | `black_holes.mass_solar` (+ `mass_class`) |
| Quasar nucleus | max(`G M / sigma^2`, radiation radius at 10 kW/m^2) | Hard | `quasars.black_hole_mass_solar`, `luminosity_w` |
| Neutron star, pulsar (magnetar if added) | galactic Hill radius between systems; in a system `max(1 T field radius, 1 g tide radius, spin-down radiation radius)` | Hard | `neutron_stars.mass_solar`, `radius_km`, `spin_period_ms`, `magnetic_field_gauss` |
| Nebula (all classes), supernova remnant | none | Pass-through with a note (chord, column, extinction) | `nebulae.density_cm3`, `extinction_av`, shape; `supernova_remnants` age, radius |
| Facility | none | Endpoint only | |

Why a sphere per kind and not a physical hazard radius: most physical radii (tidal, radiation, magnetic, accretion disk) are 1e-4 to 1e2 AU, against 1e5 AU for a star's galactic Hill sphere and 4 pc (about 8e5 AU) for a sector edge. They are below useful resolution between systems and matter only inside a system. The numbers behind each row are in section 4.

### 2.1 Below resolution, exemptions, scale hand-off

- **Galaxy and sector scale:** only galactic-scale radii matter (stars, black holes, neutron stars, large fields, large quasars). Tide, field, radiation and photon-sphere radii are under 1e-3 of a sector edge and are ignored.
- **Small spheres** (rogue planets, small fields, comets) stay obstacles when a leg crosses one, but a sphere with `r < 0.002 x leg length` is "below resolution": two tangent waypoints, no arc samples, not drawn (length change about 2e-4 percent). Inside a system the same applies below about 1e-5 of the leg.
- **Exemption by containment**, not only by "the end body": a sphere that contains the start or the goal is dropped (a planet is inside its star's galactic Hill sphere; a system can lie inside a nebula or a quasar's sphere of influence). Closed test: centre distance <= R x (1 + 1e-9). The unit tests should pin it.
- **Hand-off between scales:** the galactic path ends at the system's star (exempt sphere). The in-system leg starts where the path crosses the system's heliopause (`stars.heliosphere_radius_km`, compressed in a cloud by `compressed_heliosphere_radius`; see navigation-frames.md) and ends at the target body. The galactic Hill sphere is the exemption radius, not the hand-off radius.

## 3. Fixes to the built rule

### 3.1 Asteroid fields (NAV.51)

`keepout.py` and `queryDb.keep_out_radius` return `KeepOut(None, "none", PASS_THROUGH_NOTE)` for `asteroid_field`. Change to a hard keep-out of `radius_ly` converted to km, basis `"radius"`. The whole stored sphere is the radius, with no density-based shrinking. [nebula-and-asteroid-field-classes.md](nebula-and-asteroid-field-classes.md) section 9 suggests a 1.1 times margin as a cheap default; the planner already places waypoints at `R (1 + 1e-6)`, so no margin is needed (a tunable if Boss wants one).

This is a game rule, not physics. Expected strikes `N = n pi (r_body + r_ship)^2 L` for a 100 m ship and bodies of 1 km or more [C, `collide.py`]:

| Crossing | at belt density (4e6 km spacing) | at 1e6 km spacing |
|---|---|---|
| Belt, 1 AU | 2e-12 | 1.4e-10 |
| Field, 0.01 ly across | 3e-9 | 1.8e-7 |
| Field, 2 ly across | 3e-7 | 1.8e-5 |

A field would need to be 5.6e4 to 3.4e6 times denser than the real belt for one expected strike on a 2 ly crossing. Real belt spacing is about 1e6 km on average (4e6 km from n = 5.9e4 per AU^3 over 1 km, nebula document section 9); Pioneer 10's closest approach to any known asteroid was 8.8e6 km and Dawn's planned closest approach to a catalogued one 1.0e6 km [S]. The density label matters only for the note and a possible soft cost. Dust at speed is the one place physics supports avoidance (section 3.4 and the nebula document).

An in-system belt stays pass-through with a note: a belt lies between the inner and outer planets, so a hard keep-out would forbid most inner-to-outer legs or force a climb out of the plane, and the crossing risk is about 1e-10.

### 3.2 Supermassive black hole and quasar

Their galactic radius is 0, so their galactic Hill radius is 0 and the built floor is the event horizon: about 1.2e7 km = 0.08 AU for 4e6 solar masses. A course would graze the hole. Use the sphere of influence `G M / sigma^2`: 1.56 pc for 4e6 solar masses at 105 km/s, 1.7 pc with the generator's constant sigma = 100 km/s (which `compact_remnant.py` already stores as `system_perimeter` on the generated hole, not in the `black_holes` table). For a quasar use the larger of that and a radiation radius: luminosities in `tuning` (0.1 to 1 L_Edd, 1e38 to 1.3e40 W) at 10 kW/m^2 give 0.9 to 9 pc [C]. The region will contain many stars; any course to a star inside is exempt by containment (section 2.1) and every other course goes around.

### 3.3 Moon Hill radius (verified in the code)

Confirmed by reading the code. `Planet(...)` for a moon is created in `physics/planets.py` `generate_moons` with `planet.star` as its star and the moon's distance from its planet (AU) as `distance`. `generate_planet_properties` ends with an unconditional `update_hill_sphere(planet)`, which computes `calculate_hill_sphere(planet.distance * AU_TO_M, planet.mass, planet.star.mass)`: the moon's mass against the star's mass, at the moon's distance from its planet. `generation/validation.py` `refresh_moon_orbit` calls the same function for moons after an edit. For the Moon (384,400 km, 7.35e22 kg) this gives 888 km (smaller than its 1,737 km radius); against the Earth's mass it is 61,525 km [C, hand arithmetic from the code's formula; the project could not be imported here because `astropy` is not installed, and no database was available].

Effects beyond routing: `min_orbit_distance` for a moon is `5 x hill_radius`, and `generate_moons` adds it to the running orbit distance, so moons are spaced with the wrong number (about 4,440 km for a Moon-like moon against about 307,000 km with the right Hill radius). The outer bound (`MOON_PROGRADE_STABLE_HILL_FRACTION` of the planet's Hill radius) is correct. The facility orbit slider (`facilities.orbit_limits`) falls back to a multiple of the host radius when the sphere is below the lowest orbit, which hides the error for small moons. A moon is inside its planet's Hill sphere anyway, so routing rarely needs the value; compute it against the planet's mass (`a (m / 3 (M + m))^(1/3)`, section 8) and fix `update_hill_sphere` for moons. This is a generator bug, separate from routing.

### 3.4 Nebulae and supernova remnants

Recommendation: pass-through with a note, which is what is built, but for better reasons than the current note ("No mass is stored for it, so the course passes through it.") gives.

- **Geometry.** A nebula's drawn shape fills only 13% of its bounding sphere (6% to 19% over 30 shapes drawn with `nebula_shape.draw_shape`, seeded with `random.Random` in place of `draw.Stream`) [C, `nebfrac.py`]. A hard keep-out on the stored `radius_ly` would reserve about 7 times the real volume and force detours around empty space. If Boss ever wants a hard keep-out it should use the shape, not the sphere.
- **Gas.** Ram pressure is negligible (0.015 Pa at 0.1c in a 1e4 cm^-3 core), but kinetic power flux `0.5 rho v^3` is not: 31 W/m^2 at 0.1c in 1 cm^-3 gas, 3.1e5 W/m^2 at 0.1c in a 1e4 cm^-3 core, 3.9e7 W/m^2 at 0.5c [C; hydrogen mass 1.4 m_H]. The 50 kW/m^2 threshold is reached at 0.1c for about 1,600 cm^-3. Dust erodes about 0.5 mm of surface per 3e17 cm^-2 of hydrogen column at 0.2c (Hoang et al. 2017 [S]); crossing 100 ly of a class M core (1e4 cm^-3) is 9.5e23 cm^-2, about a million times the Alpha Centauri trip [C]. This applies only to sub-light speeds near 0.1c or more; warp and fold speeds are game physics. It is a reason to put the column density in the note, not to forbid passage.
- **Extinction.** A_V up to 100+ in dark cores blocks sensors and starlight: a navigation-difficulty note and a possible soft cost, not a collision hazard.
- **Remnants.** Hot gas (1e6 to 1e7 K) at 0.1 to 100 cm^-3, so gas hazards are small; the particle and X-ray environment of a young remnant (classes R, T, U; under about 1e4 years) is the real caution. The lethal distances in the literature (8 to 10 pc for a supernova, 91 to 200 pc for a GRB [S]) are for biospheres, not ships. Remnants have no stored shape, so the bounding sphere is the whole thing.
- **New note text.** "Passes through <name> (class X) for about N ly; column about Y cm^-2; extinction A_V about Z." The chord length comes from `NebulaShape.contains_many` sampled along the leg. A later optional "avoid dense clouds" switch (edge cost multiplied by 1 + k inside the shape, A* on a coarse grid) only if Boss asks.

## 4. Numbers behind the table

All galactic Hill radii use the code's formula (total Milky Way mass 1.15e12) at 8 kpc unless stated. [C] values are from `nav05/physics.py`.

### 4.1 Stars

| Star | Mass (solar) | Galactic Hill radius | Percent of 4 pc sector edge |
|---|---|---|---|
| M dwarf | 0.1 | 0.25 pc (0.80 ly) | 6% |
| M/K | 0.5 | 0.42 pc | 10% |
| Sun, white dwarf | 1 | 0.53 pc (1.73 ly) | 13% |
| B | 3 | 0.76 pc | 19% |
| O | 30 | 1.65 pc (5.4 ly) | 41% |
| O3 | 60 | 2.07 pc (6.8 ly) | 52% |

The galaxy's tidal (Jacobi) radius of the Sun is larger: `r_J = (G M / (4 A (A - B)))^(1/3)` with the Oort constants gives 1.39 pc, and a flat rotation curve (236 km/s at 8.2 kpc) gives 1.37 pc [C; formula and constants S, the 1.4 pc figure itself was not in the results]. The code's value is 2.6 times smaller because it uses the whole-galaxy mass (halo to the virial radius) instead of the mass inside the Sun's orbit (about 1.06e11 solar masses for a flat curve, which gives 1.2 pc as a point-mass Hill radius). [orbital-solvers-and-integrators.md](orbital-solvers-and-integrators.md) section 6.1 tabulates the point-mass Hill and Jacobi radii by galactocentric radius (1.22 and 1.37 pc for 1 solar mass at 8.128 kpc), and [galactic-potential.md](galactic-potential.md) section 10 notes that `galactic_hill_radius_km` could use the potential's enclosed mass. **Do not change this for routing.** The smaller radius is what sector placement uses to keep systems apart (`galaxy/sector.py` `hill_radius_ly`); with 1.4 pc the spheres would fill 108% of the volume at 0.1 stars/pc^3 (against 6% now) and overlap everywhere. The 6% fill is why routing around them is practical.

Other stellar options, for completeness:
- Radiation `r = sqrt(L / (4 pi F))`. At F = 10 kW/m^2: M dwarf 0.012 AU, Sun 0.37 AU, B star (1e3 solar luminosities) 12 AU, O star (1e5) 117 AU, O3 (1e6) 369 AU [C]. Parker Solar Probe's closest approach (9.86 R_sun) sees 650 kW/m^2 [C], so shielded craft go far inside this. At the F_lim = 50 kW/m^2 chosen here (a tunable): 0.17 AU for the Sun, 165 AU for an O3 star.
- Photosphere or corona floor `k R_star`; with k = 10 the Sun is 0.047 AU. The in-system star rule is the larger of the two.
- The O/B star H II region (Stromgren radius, 0.68 to 14.6 pc for Q = 1e49 photons/s at n = 1,000 to 10 cm^-3 [C]) is a nebula extent, not a ship hazard.

### 4.2 Planets and moons

- Hill radius `a (m / 3M)^(1/3)`: Mercury 2.2e5 km, Earth 1.5e6 km, Mars 1.1e6 km, Jupiter 5.3e7 km (0.36 AU), Saturn 6.5e7, Uranus 7.0e7, Neptune 1.16e8 km (0.78 AU) [C]. A hot Jupiter at 0.05 AU: 5.1e5 km (7 R_J).
- Stable satellites: Hamilton and Burns (1992, Icarus 96:43) say the stable zone scales with the Hill sphere at pericentre and retrograde orbits are stable out to slightly past the Hill radius [S]. The common numbers are about 0.49 R_H prograde and 0.93 R_H retrograde (Domingos, Winter and Yokoyama 2006) [R]; the repo uses 0.4895 (`MOON_PROGRADE_STABLE_HILL_FRACTION`). This limits where a moon can orbit; for a ship it only marks where capture and long orbits are possible.
- Laplace sphere of influence `a (m/M)^(2/5)`: Earth 9.2e5 km, Jupiter 4.8e7 km, 60% to 90% of the Hill radius [C]. No advantage over the stored value, and never mixed with it (orbital-solvers-and-integrators.md section 7).
- Fluid Roche limit `2.44 R_p (rho_p / rho_body)^(1/3)` (1.26 for rigid bodies, orbital-updates.md): 2.4 R_p for a body of the planet's own density, under 1% of a giant's Hill radius (Jupiter's is 740 R_J). A ship is rigid and strong, so Roche is a floor for parking orbits, not a hazard.
- **Default: stored Hill radius, hard keep-out.** Boss's rule, and the right order of magnitude.

### 4.3 Black holes

Schwarzschild radius `r_s = 2GM/c^2`; photon sphere 1.5 r_s; innermost stable circular orbit 3 r_s non-spinning, about 0.6 r_s prograde and 4.5 r_s retrograde at spin 0.998 [R, textbook]. Tidal acceleration across a 1 km ship equals 1 g at `r = (2 G M L / g)^(1/3)` [C].

| Class | Mass (solar) | r_s | Photon sphere | 1 g tide (1 km ship) | Galactic Hill (8 kpc) |
|---|---|---|---|---|---|
| Stellar | 10 | 30 km | 44 km | 6.5e4 km | 1.14 pc |
| Intermediate | 1e3 | 2,950 km | 4,430 km | 3.0e5 km | 5.3 pc |
| Intermediate | 1e5 | 2.95e5 km | 4.4e5 km | 1.4e6 km | 24.6 pc |
| Supermassive (Sgr A*) | 4e6 | 1.18e7 km (0.08 AU) | 1.8e7 km | 4.8e6 km | 0 (at the centre) |
| Quasar nucleus | 1e8 to 1e10 | 2 to 200 AU | | 1.4e7 to 6.5e7 km | 0 (at the centre) |

The tidal-disruption radius of a star is `R_star (M_bh / M_star)^(1/3)` (Hills 1975), equal to r_s at about 1e8 solar masses for a Sun-like star [S]; for a ship it is far smaller than the 1 g tide radius. Accretion-disk radiation at 1% of Eddington for a 10 solar-mass hole is 1.3e30 W, 21 AU at 10 kW/m^2 [C]. All of this is far below the galactic Hill and sphere-of-influence radii, so the rule for black holes is a "gravity well" rule in Boss's sense.

### 4.4 Neutron stars, pulsars, magnetars

Mass 1.1 to 2.2 solar, radius 10 to 13 km. The DB has `spin_period_ms`, `magnetic_field_gauss` (up to 1e13 G) and `pulsar_type`; there are no magnetars in the generator (1e14 to 1e15 G [S]).

- Light cylinder `c P / (2 pi)`: 67 km (1.4 ms), 760 km (16 ms), 9.5e4 km (2 s) [C; formula S].
- Dipole field `B(r) = B0 (R/r)^3` falls to 1 tesla at 260 km (1e8 G ms pulsar), 5,600 km (1e12 G), 12,000 km (1e13 G), 56,000 km (1e15 G magnetar) [C].
- 1 g tide across a 1 km ship at 1.4 solar masses: 3.4e4 km [C].
- Vacuum-dipole spin-down power `B^2 R^6 Omega^4 / (6 c^3)` is 7e27 W (ms pulsar) to 4e31 W (young, 16 ms), which at 10 kW/m^2 is 1.6 to 125 AU [C, formula R]. The beam's orientation is not stored, so `sqrt(E / (4 pi F))` is the orientation-free stand-in.
- Default: galactic Hill radius (0.59 pc) between systems and in sectors; inside a system `max(1 T field radius, 1 g tide radius, spin-down radiation radius)`.

### 4.5 Rogue planets, hypervelocity stars, comets, dust

- Rogue planet: galactic Hill radius (as built), 0.0076 pc (0.025 ly) for an Earth mass, 0.052 pc for a Jupiter mass; about 58 per local sector at the stored rate; 0.016 expected hits per 100 pc line for Earth-mass ones, below drawing resolution.
- Hypervelocity star: an ordinary star system with `runaway_speed_kms` and a velocity vector; use the star rule. It moves 0.33 ly in a century at 1,000 km/s, so ignore motion except on trips of decades.
- Comets: galactic Hill radius of a 1e12 kg interstellar comet 0.087 AU [C]; no keep-out. A coma (1e5 to 1e7 km) is a note only.

## 5. The path-bending algorithm (NAV.26)

### 5.1 Known results

- For disjoint spheres in 3D, a shortest path is straight segments tangent to spheres joined by arcs of great circles on the sphere surfaces. This is the standard picture for convex smooth obstacles [R]; the general 3D shortest-path problem among obstacles is NP-hard with an FPTAS (Kisfaludi-Bak, ISAAC 2025) [S]. Exact algorithms exist for one, two and three spheres (Chou, Hsieh and Lee, IEEE 2008; Chou 2009; a fixed-sequence n-sphere version in 2012) and an approximation scheme for disjoint spheres [S]; the papers were not read.
- Tangent lines between two spheres are continuous families in 3D, so no exact visibility graph exists; any graph needs sampled points per sphere. A lattice A* built as the reference was never better than the method below.
- **One sphere, closed form and planar.** If the straight segment A to B meets sphere (c, R), the shortest path lies in the plane through A, B and c. With `a = |A-c|`, `b = |B-c|` and `phi` the angle at c between them, the wrap exists when `acos(R/a) + acos(R/b) < phi`. Entry and exit are at angles `acos(R/a)` and `phi - acos(R/b)` from the direction of A in that plane. Length: `sqrt(a^2-R^2) + R (phi - acos(R/a) - acos(R/b)) + sqrt(b^2-R^2)`. If A, c, B are collinear any plane through the line gives the same length. The planner reproduces it to 7e-15 (compared at the waypoint radius `R (1 + 1e-6)`). Overhead for a dead-centre sphere of radius r on a leg of length L: 2e-4 at r/L = 0.01, 1.8e-3 at 0.03, 2.0e-2 at 0.1, 0.128 at 0.25, 0.54 at 0.49 [C]; about `2 (r/L)^2` for small r/L.

### 5.2 Methods compared

| Method | Result in the tests | Verdict |
|---|---|---|
| **Insert a wrap at the first hit, polish by block-coordinate descent with exact tangents** | Valid and shortest in every disjoint scene; 1 to 8 ms for up to 600 spheres | Use for NAV.26 |
| Lattice A* (K = 12 to 64 points per sphere plus tangent rings from each end) | Never shorter; equal within 0.002% in a wall of spheres, up to 3% longer with K = 16 in dense random fields; 5 ms to 1 s | Fallback or cross-check only |
| RRT*, roadmaps, exact n-sphere algorithms (Chou et al.) | Not run | Not needed: obstacles are few and convex and the local solution is closed form |

Procedure, in prototype terms (`avoid.py`; the recommended method is `method="greedy"` in `_plan_core`, reached through `plan_corridor(..., method="greedy")`; the prototype's default, `"lattice"`, is the A* reference and is not to be copied):

1. Drop spheres that contain the start or the goal. If the straight segment misses all the others, return it.
2. Loop: find the first leg that still crosses a sphere (closest-approach test; no squared large numbers), and insert a wrap using the one-sphere solution between the leg's two ends. The inserted wrap's neighbouring legs leave the sphere outward, so one wrap cannot be hit twice by the same sphere and the loop ends.
3. Polish: sweep over the wraps, recompute each wrap's tangent points from its two neighbours (exact), and accept only if the new legs and arc are clear of all spheres and the length does not rise. Drop a wrap only if the direct neighbour-to-neighbour segment is clear of all spheres. Merge adjacent wraps on the same sphere. Repeat until the length stops falling.
4. A first version let the polish drop and re-insert wraps freely; it oscillated and left 5% to 20% of dense scenes invalid. The validated polish above made every disjoint test valid.
5. Broad phase for the full data set: only spheres whose centre is within R + W of the straight segment go to the planner, W = max(5% of the leg, 3 x the biggest radius in that tube). After planning, check the result against all spheres and double W on any miss (at most 6 rounds). In the DB this is NAV.25's corridor query.

### 5.3 Key function bodies

The hit test uses the closest-approach form. The discriminant form `B^2 - 4A(C - R^2)` of "Orbital Update Collision Detection.md" loses digits when R is far smaller than the coordinates (the same cancellation [collisions-and-mergers.md](collisions-and-mergers.md) section 2.1 documents for that code). For a ship leg against a static sphere it is the same quadratic with `Delta v = 0` for the sphere.

```python
def _closest(P, Q, C):
    d = Q - P; a = float(d @ d); f = C - P
    tu = (f @ d) / a                          # unclipped closest-approach parameter
    t = np.clip(tu, 0.0, 1.0)
    dist = np.linalg.norm(f - np.outer(t, d), axis=1)
    return t, dist, tu, math.sqrt(a)

def seg_hits(P, Q, C, R):                     # open ball: tangent grazing is not a hit
    if float((Q - P) @ (Q - P)) == 0.0:
        return np.zeros(len(R), bool)
    t, dist, tu, dl = _closest(P, Q, C)
    return dist < R

def seg_first_hit(P, Q, C, R):                # sphere entered first along P->Q, or -1
    if float((Q - P) @ (Q - P)) == 0.0:
        return -1
    t, dist, tu, dl = _closest(P, Q, C)
    hit = dist < R
    if not hit.any():
        return -1
    half = np.sqrt(np.maximum(R * R - dist * dist, 0.0)) / dl
    t_in = np.where(hit, np.maximum(tu - half, 0.0), np.inf)
    return int(np.argmin(t_in))

def tangent_wrap(a, b, c, R):
    """Exact shortest a -> b around ONE sphere. None if the segment clears it,
    else (entry, exit, arc_angle). Lives in the plane through a, b and c."""
    u, v = a - c, b - c
    da, db = np.linalg.norm(u), np.linalg.norm(v)
    if da < R * (1 - 1e-9) or db < R * (1 - 1e-9):
        return None
    e1 = u / da
    w = v - (v @ e1) * e1; wn = np.linalg.norm(w)
    if wn < 1e-12 * db:                       # collinear: any perpendicular will do
        k = int(np.argmin(np.abs(e1))); t = np.zeros(3); t[k] = 1.0
        w = t - (t @ e1) * e1; wn = np.linalg.norm(w)
    e2 = w / wn
    phi = math.atan2(v @ e2, v @ e1)          # in (0, pi]
    alpha_a = math.acos(min(1.0, R / da)); alpha_b = math.acos(min(1.0, R / db))
    if alpha_a + alpha_b >= phi:
        return None
    ang_e, ang_x = alpha_a, phi - alpha_b
    e = c + R * (math.cos(ang_e) * e1 + math.sin(ang_e) * e2)
    x = c + R * (math.cos(ang_x) * e1 + math.sin(ang_x) * e2)
    return e, x, ang_x - ang_e

def repair(items, C, R, max_inserts=200):     # items: [("pt", p) | ("wrap", i, entry, exit)]
    Rt = R * (1 - 1e-7)                       # hit tests use R(1-1e-7), waypoints sit at R(1+1e-6)
    for _ in range(max_inserts):
        for k in range(len(items) - 1):
            a = items[k][1] if items[k][0] == "pt" else items[k][3]
            b = items[k + 1][1] if items[k + 1][0] == "pt" else items[k + 1][2]
            j = seg_first_hit(a, b, C, Rt)
            if j < 0: continue
            res = tangent_wrap(a, b, C[j], R[j] * (1 + 1e-6))
            if res is None: continue
            items.insert(k + 1, ("wrap", j, res[0], res[1]))
            break
        else:
            return items                      # no leg crosses a sphere
    return items

def polish(items, C, R, sweeps=60):           # block-coordinate descent, exact tangents
    for _ in range(sweeps):
        before = chain_length(items, C, R)
        for k in range(1, len(items) - 1):
            if items[k][0] != "wrap": continue
            i = items[k][1]
            a, b = end_of(items[k - 1]), start_of(items[k + 1])
            res = tangent_wrap(a, b, C[i], R[i] * (1 + 1e-6))
            if res is None:                   # a-b clears this sphere
                if seg_first_hit(a, b, C, R * (1 - 1e-7)) < 0:
                    drop_wrap(items, k)       # only if a-b clears every sphere
                continue
            cand = ("wrap", i, res[0], res[1])
            if wrap_ok(a, cand, b, C, R) and new_length(cand) <= old_length + 1e-12:
                items[k] = cand               # clear of every sphere and never longer
        merge_adjacent_wraps_on_same_sphere(items)
        if before - chain_length(items, C, R) < 1e-12 * before: break
    return items

def repair_polish(items, C, R, sweeps=60, rounds=40):
    for _ in range(rounds):
        items = polish(items, C, R, sweeps)
        new = repair(items, C, R)
        if len(new) == len(items): return items
        items = new
    return polish(items, C, R, sweeps)
```

`end_of`, `start_of`, `drop_wrap`, `wrap_ok` (both legs and a 0.05 rad sampling of the arc clear of every other sphere), `chain_length` (legs plus `radius x arc angle`), `new_length`/`old_length` (that wrap with its two neighbours, before and after) and the merge helper are small functions in `avoid.py`; `plan_corridor` wraps the whole thing in the broad phase of step 5. The prototype needs only `math` and numpy, which the project already depends on; a pure-Python version was not tested.

### 5.4 Measured results

All [C]. Scenes are random fields of spheres in a 100-unit cube unless stated; "overhead" is path length over straight line minus 1; "valid" means a dense-sample check found no point inside any sphere by more than 1e-7 relative.

| Scene | Cases | Valid | Overhead mean / p95 / max | Time mean / p95 |
|---|---|---|---|---|
| Disjoint, 5% fill, r 2 to 10 | 40 | 100% | 0.17% / 0.36% / 3.9% | 1.4 / 10 ms |
| Disjoint, 12% fill, r 3 to 15 | 30 | 100% | 0.13% / 0.91% / 1.6% | 2.8 / 16 ms |
| Disjoint, 25% fill, r 4 to 16 | 20 | 100% | 0.60% / 3.0% / 3.0% | 8 / 21 ms |
| Two-layer wall of spheres with 2-in-10 gaps | 12 | 100% | 0.09% / n/a / 0.38% | 3 ms (lattice 190 ms, same length within 0.002%) |
| Star field (0.114/pc^3, 7% fill, 24,568 spheres): hop 10 pc | 60 | 100% | 0.20% / 1.2% / 2.2% | 19 / 49 ms |
| same, hop 30 pc | 60 | 100% | 0.10% / 0.45% / 0.64% | 40 / 82 ms |
| same, hop 90 pc | 20 | 100% | 0.07% / 0.25% / 0.34% | 151 / 273 ms |
| Overlapping spheres, 25% fill | 25 to 40 | 84% to 90% (lattice no better: 84%) | 1.4% mean | 4 ms |
| Overlapping, 50% fill | 40 | 72% | | 3 ms |

- The star-field times are dominated by whole-array distance tests over 24,568 spheres with no spatial index. With NAV.25's corridor query (a few hundred spheres) expect 1 to 10 ms.
- Hops of 10 pc bend in 47% of cases and hops over 30 pc in 88% or more: about one wrap per 15 pc of line (0.7, 1.8, 4.1 and 6.1 wraps at 10, 30, 60 and 90 pc). Overhead stays under about 0.2% on average, so adding keep-out does not change which stops are chosen, and travel times grow by a few tenths of a percent at most. Ignoring a sphere of radius r on a leg of length L would change the length by about `2 (r/L)^2` (under 0.02% for r < 0.01 L), but small spheres are still stepped around so the path is correct.
- Robustness: a chain of 40 tiny spheres on the line (radius 1e-6 to 1e-2 of a 100 pc leg, or 1e-4 to 1e-3 pc on a 2e4 pc leg) gave the same result to 1e-15 relative with the whole scene shifted by up to 2e4 pc. Work in parsecs (galaxy), light-years (sector) or km (system) with doubles.
- Ten edge cases pass (`t19.py`): start equal to end; start inside a sphere; end on a sphere surface; start inside a sphere that overlaps another on the way; collinear start-centre-end; two tangent spheres with the line through the contact point (goes straight); two spheres 0.01 apart (goes through the gap); a sphere near but not touching the line; both ends in one big sphere; a zero-radius sphere.

### 5.5 Degenerate and hard cases

- **Start or end inside a sphere:** drop that sphere (exempt). **Both ends inside the same sphere:** straight line.
- **Overlapping spheres:** the arc of one sphere may enter another; validity fell to 72% to 90% at heavy random overlap. Real data should rarely overlap because system Hill spheres are placed disjoint (`galaxy/sector.py`), but black holes, rogue planets and fields are not covered by that rule. Handling: detect clusters (union of overlapping pairs); if the plan still crosses a sphere after the loop, replace the cluster by the smallest enclosing sphere (`merge_overlaps`, exact for pairs) and replan; if that sphere swallows an end (it did in every test at 25% fill, because random clusters percolate), return the straight line with the warning `overlap_fallback` and not a wrong path. An overlap census at generation time or in a test shows how often this fires.
- **Many spheres:** the loop is linear in wraps; broad phase by corridor; cap the width doublings (6) and the number of wraps (a few hundred) and return a flagged straight line if exceeded.
- **Chosen side:** the method keeps the first wrap's plane. It never lost to the lattice search in random scenes and was within 0.002% in the two-layer wall, but it is not proven globally optimal; the UI should say "short detour", not "shortest".
- **Rounding:** waypoints sit at `R (1 + 1e-6)`, hit tests use `R (1 - 1e-7)`, so tangent legs never fail on rounding.

## 6. Moving bodies inside a system (NAV.27)

Real planet data (a, mass, circular coplanar orbits) from Mercury to Neptune, random phases, legs Earth to Neptune, Mercury to Uranus and Venus to Saturn [C, `moving.py`, `moving2.py`; 300 and 200 draws]:

| Ship speed | Static t0 positions enter a Hill sphere (more than 1% of R_H) | Iterated (obstacle at time of pass) |
|---|---|---|
| 1c, 10c, 100c (warp 1 and up) | 0% (max 0.00 R_H), also 0% at 0.1c | 0% |
| 0.01c | 0.5% to 5.5% (max 0.10 R_H) | 0% |
| 0.003c | 1.5% to 2.5% (max 0.26 R_H) | 0% |
| 0.001c (300 km/s) | 1.5% to 2.5% (max 0.97 R_H) | 0% |
| 0.0003c (90 km/s) | 0.5% | 0% to 1% (does not converge near planet speeds) |

Ship speed at warp 1 is 3e5 km/s against 5 to 48 km/s for planets, and outer Hill radii are 5e7 to 1.2e8 km. Inner planets are the only close call: Mercury's Hill radius is 2.2e5 km, and a ship at 1c needs about 140 s to cross 0.28 AU, during which Mercury moves 6,600 km (3% of R_H).

Recommended procedure (no Lambert solution, no ephemeris iteration):

1. Plan with every body at the departure time t0 (`body_positions.positions_at`); enough for warp or fold.
2. The destination body moves while the ship flies: aim at `T(t_arrival)` with `t_arrival = t0 + L/v`, a fixed point that converges in 2 to 3 steps because `u/v` is below 2e-4.
3. If a speed below about 0.01c is ever offered, do 3 to 6 passes: set each obstacle's time to `t0 + s_i / v` (arclength along the previous path to its closest approach) and replan.
4. A comet on a parabolic or hyperbolic orbit is handled the same way (`positions_at` already does it).
5. Default departure time: now (NAV.27 text).

At galactic scale nothing moves enough at warp (a star at 200 km/s drifts 0.067 ly a century; a hypervelocity star at 1,000 km/s, 0.33 ly), so NAV.27 stays system-only.

Gravity itself does not bend the ship in this model: warp and fold speeds are game physics and the keep-out is a rule. The flyby formula of "Computational Astrodynamics.md" (turning angle `delta = 2 arcsin(1 / sqrt(1 + (b v^2 / mu)^2))`, correct) shows why that is harmless for routing. Passing the Sun (mu = 1.327e11 km^3/s^2) at 0.17 AU: 0.96 rad at 100 km/s, 1.2e-3 rad at 3,000 km/s, 1.2e-5 rad at 0.1c, 1e-7 rad at 1c [C]. Gravitational deflection matters only for a ship at planetary speeds, which no speed in the warp and fold tables is.

## 7. How course lines are exposed (NAV.28, NAV.5, NAV.17, NAV.21)

Add `adjusted` to the NAV result and to `galaxy_course` beside the existing `points` and `direct`. Positions use the frame already in use (galaxy parsecs for the Galaxy Map, sector-local light-years, system km), declared in `unit` and `frame`.

```json
"adjusted": {
  "frame": "galactic", "unit": "pc",
  "length": 40.2966, "straight_length": 40.1123, "extra": 0.1842,
  "waypoints": [
    {"p": [0,0,0],               "kind": "start",       "ref": null},
    {"p": [11.848,-0.55,-0.755], "kind": "tangent_in",  "ref": "star:101"},
    {"p": [12.222,-0.54,-0.757], "kind": "tangent_out", "ref": "star:101"},
    {"p": [40,3,0],              "kind": "end",         "ref": null}
  ],
  "legs": [
    {"type": "line", "from": 0, "to": 1},
    {"type": "arc",  "from": 1, "to": 2, "center": [12,1,0.5], "radius": 2.0, "obstacle": "star:101"},
    {"type": "line", "from": 2, "to": 3}
  ],
  "polyline": [[0,0,0], "..."],
  "avoided": ["star:101"],
  "passed_through": [{"ref": "nebula:55", "note": "...", "inside_length": 12.4, "column_cm2": 3.1e21}],
  "ignored_below_resolution": ["rogue_planet:9001"],
  "ruleset": {"keepout": 1, "planner": 1},
  "warnings": []
}
```

| Field | Meaning |
|---|---|
| `frame`, `unit` | `galactic` + `pc`, `sector` + `ly`, or `system` + `km`; every position in the object is in this frame |
| `length`, `straight_length`, `extra` | Adjusted path length, the straight line, and the difference; same unit. NAV.2's tables use `length / speed` unchanged |
| `waypoints` | `{p, kind, ref}`; `kind` is `start`, `tangent_in`, `tangent_out` or `end`; `ref` is the object reference (`star:101`) of the sphere a tangent point sits on |
| `legs` | The exact form: a `line` between two waypoint indices, or an `arc` with `center`, `radius` (the radius used, `R (1 + 1e-6)`) and `obstacle`; an arc is the shorter great-circle arc between its two waypoints |
| `polyline` | What the map draws; arcs sampled so the chord sagitta is at most 1% of the sphere radius (about 12 segments for a half circle; 1 to 2 for a below-resolution sphere); vertices sit on or just outside the sphere |
| `avoided` | References of the spheres the path bends around |
| `passed_through` | Pass-through objects crossed, with the note text, length inside and column density (section 3.4) |
| `ignored_below_resolution` | Spheres smaller than 0.002 of the leg (section 2.1), listed but not drawn |
| `ruleset` | Versions of the keep-out table and the planner, so an old saved course can be told from a changed rule |
| `warnings` | Machine-readable codes: `overlap_fallback`, `wrap_limit`, `unknown_space` (an ungenerated sector has no objects to avoid; the NAV.12 flag already marks it) |

- The straight line stays as `direct`; the system-to-system hop chain as `route`. Per-hop adjusted forms go inside the route entries and the route's total is their sum. NAV.23's switch (Direct, Route, Adjusted) is three keys.
- The course readout (bearing and mark, [navigation-frames.md](navigation-frames.md)) is worked out per leg in that leg's frame; for an adjusted course the heading to leave on is the first leg's, not the direct line's.
- Saved courses (NAV.17): store `waypoints`, `legs`, the `ruleset` version and, for each obstacle, `ref` and the radius used (a few hundred bytes per wrap). Opening a saved course recomputes it and says what changed: an object moved (the correlative update), a regenerated system vanished, or the keep-out rule changed. Without the radius used there is no way to tell which.
- NAV.21 (fit the view): take the bounding box of `polyline`, not of `waypoints`, since an arc can bulge past its tangent points.
- MAP.121: course geometry is in the galaxy frame, independent of how neighbours are dimmed or faded; draw the course on a layer that is not faded with the sectors in front of the camera.
- Browser twin: `sectors_along_segment` has one; the planner can stay server-side (the result is data), so no JavaScript port is needed unless Boss wants live re-planning while dragging endpoints.

## 8. What carries over from the uploaded documents, and what does not

- **"Computational Astrodynamics.md", Hill sphere.** The formula `d (m / 3M)^(1/3)` is what `orbits.calculate_hill_sphere` implements and what the stored `hill_radius_km` uses. Against the exact L1 distance of the restricted three-body problem it is +0.3% for Earth and +2.4% for Jupiter, but +6% for the Moon against the Earth, +14% at m/M = 0.1 and **+39% for equal masses**; with `(m / (3 (M + m)))^(1/3)` the equal-mass error is +10% ([orbital-solvers-and-integrators.md](orbital-solvers-and-integrators.md) section 1, row 5). For keep-out this matters only for moons (use `M + m`, section 3.3) and a binary pair's mutual sphere; planets keep the stored value. The document's `d` is the instantaneous separation, while the stable sphere is set at pericentre `a (1 - e)`; the stored value uses `a`, the generous side for a ship, so it stays.
- **Same document, Roche limit and hyperbolic deflection.** `d_roche = 2.44 R1 (rho1/rho2)^(1/3)` is the fluid limit and holds for a fluid secondary; it is not a ship hazard (section 4.2). The deflection formulas are correct as written and are used in section 6.
- **Same document, sector frame and galactic potential.** It describes 4 pc cubic cells with Morton keys; planetGen uses 4 pc ring, layer and slot cells (orbital-solvers-and-integrators.md row 12), so courses use the existing frames: galaxy parsecs, sector-local light-years, system km. Its galactic potential gives an enclosed mass for a Jacobi radius; the code deliberately uses the total mass (section 4.1), and [galactic-potential.md](galactic-potential.md) records the option to change it. This design does not.
- **"Orbital Update Collision Detection.md".** Its swept-sphere quadratic `S(tau) = A tau^2 + B tau + C` is the mathematics of a leg against a static sphere. The planner uses the closest-approach form for the cancellation reason in section 5.3; straight-line sweeps are exact for ship legs, and the curved-path limitation (orbital-solvers-and-integrators.md section 7) applies only to moving bodies, which NAV.27 evaluates at the time of pass. Its other defects, and the merger radius, are in [collisions-and-mergers.md](collisions-and-mergers.md) section 2.

## Evidence notes

[S] seen in a search result (snippet text): the Hamilton and Burns stability statement; the Jacobi radius formula with Oort constants (A about 14.8, B about -12 km/s/kpc); belt spacing about 1e6 km, Pioneer 10 closest 8.8e6 km, Dawn closest catalogued 1.0e6 km; the tidal-disruption radius and Hills mass (about 1e8 solar masses); Hoang et al. dust erosion column; magnetar fields and the light-cylinder formula; the existence of exact three-sphere and approximate multi-sphere shortest-path algorithms and the NP-hardness statement for 3D obstacles.

[C] computed: all galactic Hill radii, Jacobi comparisons, filling fraction and hit-rate figures, black hole radii, tidal and radiation radii, neutron star radii and spin-down power, Stromgren radii, belt strike table, kinetic power flux table, nebula shape volume fraction (13%), all planner results and timings, the moving-bodies tests, the robustness and edge-case tests, the flyby deflection figures (section 6) and the moon Hill-radius figures (hand arithmetic from the code's formula). Not computed against a real database: no generated galaxy was available (importing the generator needs `astropy` and `nltk`, which were not installed), so densities and masses come from `tuning.py` and per-sector counts are model-level, not measured.

[R] recalled, to verify when paper or arXiv access is allowed:
- The 0.4895 and 0.9309 Hill fractions attributed to Domingos et al. 2006 (a result listed both numbers with the labels swapped, so the pairing needs a check).
- Photon sphere 1.5 r_s and the innermost-stable-orbit values (3 r_s, spin dependent).
- The vacuum-dipole spin-down formula.
- The main-belt count of about 1e6 bodies over 1 km and the 18 AU^3 torus estimate.
- That shortest paths around spheres are tangent segments plus great-circle arcs (standard for convex obstacles, but not confirmed by a search result).
- The M-sigma relation (not used; the generator's constant sigma = 100 km/s is).

## Sources

Search results, not fetched pages:
- Hamilton and Burns 1992, Icarus 96:43, reprint: https://pages.astro.umd.edu/~dphamil/research/DPHreprints/HamBurns92.pdf (snippet only)
- Stable satellites around giant planets (Domingos et al. 2006): https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-2966.2006.11104.x (title only; the 0.4895 and 0.9309 values were not confirmed)
- Jacobi radius with Oort constants: https://arxiv.org/pdf/1801.04278, https://arxiv.org/pdf/1904.05896, https://en.wikipedia.org/wiki/Oort_constants
- Asteroid belt spacing, Pioneer 10 and Dawn: https://www.astronomy.com/science/how-do-spacecraft-safely-navigate-the-asteroid-belt/, https://science.nasa.gov/mission/dawn/faq/, https://lucy.swri.edu/MainBeltDensity.html, https://en.wikipedia.org/wiki/Pioneer_10, https://www.scientificamerican.com/article/in-science-fiction-movies/
- Tidal disruption radius and Hills mass: https://link.springer.com/article/10.1007/s11214-021-00818-7, https://pubs.aip.org/physicstoday/article/67/5/37/414728/The-tidal-disruption-of-stars-by-supermassive, https://academic.oup.com/mnras/article/455/1/859/984391
- Supernova remnant radii and phases: https://arxiv.org/pdf/1701.05942, https://pages.astro.umd.edu/~rmushotz/ASTR480/A480_supernova_remnants_2016_lec3.pdf, https://arxiv.org/pdf/1604.04395
- Relativistic craft and the interstellar medium: https://arxiv.org/pdf/1608.05284 (Hoang et al. 2017, ApJ 837:5), https://arxiv.org/abs/2307.12160
- Magnetar fields and the light cylinder: https://solomon.as.utexas.edu/magnetar.html, https://arxiv.org/pdf/1804.05343
- Supernova and GRB lethal distances: https://iopscience.iop.org/article/10.1086/346127 (Gehrels et al. 2003), https://arxiv.org/pdf/1702.04365
- Shortest paths among spheres: https://ieeexplore.ieee.org/document/4811457/, https://www.researchgate.net/publication/255640870_A_Base_Algorithm_for_Computing_the_Shortest_Path_of_Spheres, https://cerv.aut.ac.nz/wp-content/uploads/2015/08/MItech-TR-75.pdf, https://drops.dagstuhl.de/storage/00lipics/lipics-vol359-isaac2025/LIPIcs.ISAAC.2025.46/LIPIcs.ISAAC.2025.46.pdf

Repo files read: `docs/TODO.md` (NAV entries), `course-routing.md`, `navigation-frames.md`, `nebula-and-asteroid-field-classes.md`, `orbital-solvers-and-integrators.md`, `galactic-potential.md`, `collisions-and-mergers.md`, "Computational Astrodynamics.md", "Orbital Update Collision Detection.md", `src/planetgen/galaxy/keepout.py`, `db/query.py` (`keep_out_radius`), `physics/planets.py`, `physics/orbits.py`, `physics/body_positions.py`, `generation/star.py`, `generation/binary.py`, `generation/validation.py`, `generation/phenomena/compact_remnant.py`, `galaxy/nebula_shape.py`, `population/facilities.py`, `tuning.py`, `db/schema.sql`, `web/nav_page.py`.
