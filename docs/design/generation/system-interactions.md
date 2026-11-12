# Everything in a system affects everything else in the system

The design principle behind every physics item of the generator, and a register of the code paths that still count only the host star or only one pair of bodies.

Informs: GEN.1441, GEN.1442, GEN.1443, GEN.1444, GEN.1445, GEN.1446, GEN.1447, GEN.1448, GEN.1449, GEN.1450, GEN.1451, GEN.1452, GEN.1453, GEN.1454, GEN.1455, GEN.1456, GEN.1457, GEN.1458, GEN.1459, GEN.1460, GEN.1461, GEN.1462

Status: plan, 2026-11-08; the principle is Boss's, the register is a first pass by reading the code (a search for host-only and pair-only effects), not a proof that nothing else is missing

## The principle

Boss, 2:58 PM EST on 2026-11-08: "Everything in a system effects everything else in the system."

Every body of a system (stars and other star-like objects, planets, moons, belts and asteroids, comets, rings) acts on every other body of the same system through its 3D position, averaged over the orbits: light, UV, X-ray, flares and particle dose; gravity as tides and as secular, resonant and scattering perturbations; wind and magnetic fields; dust and impacts. An effect that reads the host star alone, or one pair of bodies alone, is a simplification. It may stay only where it is stated and registered in GEN.725 with a tier, a tweak level and a reason (a test shows the missing terms are below a stated size).

## Rules for new items

1. Say which sources and which receivers the effect covers. The default is all bodies of the system on all bodies of the system.
2. Use the 3D positions and the system bounds already stored, averaged over the orbits (mean of 1/d^2 for light, mean of 1/d^3 for tides); do not use the host's distance as a stand-in for the distance to every source.
3. A lone star with no other body must give the numbers of before the change; a test says so.
4. Where a term is dropped (shadowing, higher-order resonances, an unknown source), say so in the item and register the simplification in GEN.725.
5. Every such change ships with a before and after census that shows how many bodies move (GEN.1445 for radiation, GEN.1448 for tides and secular effects).
6. Every effect fills the effect ledger (GEN.1465): for each receiving body, it calls `system.effect_ledger.add(effect, receiver, source, value)` once for every source it computes, however small the value (the ledger stores the share of the largest source, `physics/effect_ledger.py`). A source that is dropped from the effect for a registered reason is still noted in the ledger, so the census can show when it starts to matter. Light and UV keep the sources over a share of 1e-3 (0.1 percent) as rows; the five gravity effects keep only the largest source and a summary and are recomputed on demand from stored masses and orbits (point masses, a few milliseconds a system). A number that came from a random draw or a past run and cannot be recomputed (spin-down state used, lock time, obliquity before drift, eccentricity before and after perturbation, Kozai cycle count and start age, the formula version) is stored on the body by the item that builds that effect (`reports/gen1465-ledger-store-vs-compute-2026-11-09.md` section 4).

## Register of host-only and pair-only code paths (first pass)

| Effect | Where it is today | Item |
|---|---|---|
| Mean light, temperature, climate, albedo, habitable zone | Host star alone; only the size of the light swing counts the other stars | GEN.1441 |
| UV, X-ray, flare and particle dose, ozone, atmosphere loss | Host star alone | GEN.1442 |
| Sky and day text | Host star alone | GEN.1443 |
| Light of the dark members of the orbit tree (brown dwarf, neutron star, black hole with a disk fed by a companion's wind) | Built (`physics/dark_light.py`, rows DL1 to DL7 of section 6.259): they are sources of the radiation field like any star, in the temperature, the dose and the sky text. They give bolometric light only (no XUV, flare or particle terms) and are not stored in the effect ledger yet, because they have no stored row id | GEN.1501 |
| Belt, asteroid and comet temperature, ice survival and dose | Built for temperature, taxa, the size law (icy or rocky) and comet activity: the ring's or the perihelion's light counts every star (GEN.1444, rows RF6 to RF8), and a comet's activity distance scales with the host's luminosity as `sqrt(L / Lsun)` (GEN.1491, rows CA1 to CA3). Belts and comets store no surface dose or temperature, so there is no stored dose to change; a dose for them would be a new stored field | GEN.1444 |
| Tide raisers of a planet (ocean tide, spin, locking) | Ocean tide built (GEN.1446 part 1, `physics/tide_sources.py`): the host and companion stars, the sibling planets, the own moons or planet and the sibling moons, by the mean of `1 / d^3`, each noted in the effect ledger (stand-ins TS1 to TS5). Locking and spin-down still count the host and the own moons only | GEN.1446 |
| Tidal heating and tidal locking of planets and moons | The host and the body's own moons; the eccentricity heating reads the host mass alone | GEN.1446 |
| Kozai and secular pumping | Used for stellar triples only (physics/stability.py), not for planets or moons | GEN.1447 |
| Moon stability, Hill limit and moon spin-down | The planet's Hill sphere from the star alone; torques from the star and the planet's moons only | GEN.1449 |
| Planet chain stability and secular exchange | Mutual Hill spacing of neighbours and, since GEN.1450, the Laplace-Lagrange exchange of eccentricity and inclination among all the planets of a star (`physics/secular_chain.py`, kept inside the spacing); the companion star is in the spacing at its 3D place (GEN.1460) and its secular pumping is the Kozai screen's (GEN.1447), which reads `precession_rates`; stand-ins SC1 to SC5 | GEN.1450 |
| Belt and disk edges | Built: a companion star cuts the disk and its belts at its truncation radius before their mass is worked out (`physics/disk_truncation.py`, rows DT1 to DT4). Every planet's chaotic zone is already inside the keep-out (DT3); a companion's secular excitation of a belt is not built (DT4) | GEN.1451 |
| Stellar wind, magnetopause, radiation belts | The host star's wind alone | GEN.1452 |
| Impact and comet flux | The host star's debris disk and its own belts; no scattering by other giants, no companion-stirred cloud | GEN.1453 |
| Obliquity and spin-orbit secular drift | Torques of the host star and own moons; no sibling planet or companion | GEN.1454 |
| Trojan and Lagrange stability | One planet with its star | GEN.1455 |
| Ejection and scattering | The scatterer of a thrown or swept planet is the planet of the highest Safronov number among all whose reach holds it, not the recorded neighbour alone; a companion star ejects nothing, since planets form inside its stability limit and the eccentric-giant kick tail is not modelled (rogue-ejection-kick.md) | GEN.1456 (built) |
| Comet orbits | Two-body about the host; a Jupiter-family comet of a star with a giant has a Tisserand parameter about its most massive giant between 2 and 3 (`physics/comet_control.py`), drawn again until it does, else it is of another family; the companion's tide is the ceiling cut already in place and the galactic tide on the far cloud is GEN.1466 | GEN.1457 (built) |
| Seasonal swing, eclipses and a moon's extra night | One sun's declination; a moon darkened by its planet only | GEN.1458, GEN.1461, GEN.1462 |
| Settled eccentricity | Tidal wearing-down by the host; no pumping set against it | GEN.1459 |
| Orbit spacing and Hill checks | Planets: held to every planet and belt inside them by mutual Hill radii at apoapsis and periapsis (clearance ledger, GEN.1014, GEN.1023). Moons: every moon inside (`moon_clearance_au`, GEN.1460). Wide pairs: every planet of one star against every planet of the other (`cross_star_violation`, GEN.1460). The circular Hill radius of `min_orbit_distance` stays, registered as SP1 (the conservative one); the companion's tide on a moon's limit is GEN.1449, the secular pumping of an outer planet GEN.1447 and GEN.1450 | GEN.1460 (done), GEN.1449, GEN.1450 |

The full audit, with 35 entries, files and guessed sizes, is /mnt/project-files/reports/host-only-audit-2026-11-08.md (Lane 05). It lists what it did not read: planet mass and giant draw from the host's disk (a companion's truncation of the disk), hydrosphere and interior heat, web page text, and companion placement; a second pass is needed there.
