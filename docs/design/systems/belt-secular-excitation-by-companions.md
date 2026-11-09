# A companion star's secular excitation of a belt (GEN.1492, design)

Status: design for review by TODO; nothing is built. Written by Lane 22 on 2026-11-09. GEN.1451 (the companion's cut of a disk and its belts, `physics/disk_truncation.py`) and GEN.1450 (the planets' secular exchange, `physics/secular_chain.py`) are on main.

## The question

A belt's members get their eccentricity and inclination from one Rayleigh draw whose width depends only on how far out the belt is (`belt_roster.rms_eccentricity`, `E_OVER_I`), whatever stars the system holds. A belt a few AU from a star with a companion of a few hundred AU is pushed by the companion's tide: over many orbits the companion's pull averages into a slow exchange that raises every member's eccentricity to a **forced** value (the secular forcing of Heppenheimer 1978) that depends on where the belt stands in the companion's orbit, and that the free, drawn eccentricity then adds to as a vector. Today a belt in a close binary looks the same as one round a lone star. GEN.1451 left this out on purpose (row DT4 of section 6.246).

## What the code does today (read, not assumed)

- `belt_roster.draw_roster` draws each stored member's `e`, `i`, node, periapsis and anomaly from the belt's own stream: `e = rayleigh(u_e, e_rms)`, `i = rayleigh(u_i, e_rms / E_OVER_I)`, node and periapsis uniform. A scatter population (`belt_populations.orbit`) may fix the `e` or scale it. `Belt._draw_members` stores them in the member rows (GEN.1070). Nothing else in the generator reads `e_rms`.
- `radiation_field.companions_of(star)` gives each other star as a `Companion` (`axis_au`, `eccentricity`, the star with its mass). A subtree of stars is one point; the orbits are coplanar (the module says so). No orbit of a star carries a mutual inclination.
- `secular_chain.matrices(a_au, masses_kg, host_mass_kg)` gives the apsidal matrix `A` of a star's planets from Laplace coefficients, and `precession_rates` their apsidal rates with relativity. It does not include a companion star; the Kozai screen (`stability.kozai_verdict`, GEN.1447) compares a planet's rates against the companion's.
- `disk_truncation.star_disk_limit_au` already ends a belt where the companion cuts the disk or no orbit holds, so a belt that exists is well inside the companion's orbit (`alpha = a / a_c` below about 0.3).

## The physics

For a test body at `a` round a host of mass `M` with a companion of mass `m'` at `a_c`, eccentricity `e'`, the first-order secular theory (Murray and Dermott, *Solar System Dynamics* sections 7.2 and 7.3 [R3220]; Heppenheimer 1978) gives, with `alpha = a / a_c`, `n` the body's mean motion:

    A_c = (n / 4) (m' / M) alpha^2 b_{3/2}^{(1)}(alpha)  ~  (3/4) n (m'/M) alpha^3 (1 - e'^2)^(-3/2)
    B_c = (n / 4) (m' / M) alpha^2 b_{3/2}^{(2)}(alpha) e' ~  (15/16) n (m'/M) alpha^4 e' (1 - e'^2)^(-5/2)

`A_c` is how fast the companion turns the body's periapsis round, `B_c` the forcing along the companion's periapsis. The body's eccentricity vector is a **free** part, of length `e_free`, circling at `g = A_total` round the **forced** part of fixed length

    e_forced = B_c / A_total  ~  (5/4) alpha e' / (1 - e'^2)      for A_total = A_c alone,

and the planets of the star, general relativity and other companions add their own apsidal rates to `A_total` (`secular_chain.precession_rates` for the planets, as a body of no mass among them), which **dilutes** the forcing: a belt beside a heavy giant turns its periapsis round faster than the far companion can force it and keeps its own eccentricity. Several companions force independently and their forced vectors add (linear theory); a subtree is one point.

The member's eccentricity at the star's age is `e = | e_free exp(i varpi_free) + e_forced exp(i varpi_c) |`, where `varpi_free` is the member's drawn periapsis and `varpi_c` the periapsis direction of the companion's orbit, one angle for the system. The semi-major axes are unchanged.

## What it leaves out (stated)

- Inclination. The orbits of stars are coplanar in the generator, so the forced inclination is zero; a belt tilted by a companion needs a mutual inclination of the star orbits, which does not exist (propose a separate item). Kozai-Lidov pumping needs the same and stays with GEN.1447.
- First order in `e'` and `alpha`: held to `alpha < 0.3` (which the cut ensures) and `e_forced < 0.3`; beyond it the member keeps its drawn orbit.
- The orbits of the companion are as they stand today, not as they were while the disk lasted; resonances with the companion and the evection are not modelled.
- Planets are not forced by the companion here (GEN.1450 has no stellar term): that is a separate item and would use the same `A_c`, `B_c`.

## Build in slices (each small, each its own PR)

0. **`physics/belt_forcing.py`, calibration (no output change).** `forced_eccentricity(a_au, host_mass_kg, companion, planet_rates)` from exact Laplace coefficients (`resonant_chains._laplace`), `companion_rate`, and `combined_orbit(e_free, varpi_free, forced)`. A census over generated systems with companions: the share of belts whose forced eccentricity exceeds the drawn rms, the dilution by giants, and how many members the combination would move past the stability limit. Stand-in rows: FE1 the first-order law (L0, recalled), FE2 the companion's periapsis direction drawn per system (L1), FE3 inclination taken as zero (L0), FE4 the validity limits. Nothing is wired.
1. **Wire it into `Belt._draw_members` (lock item).** The members of a belt of a star with companions get the combined `e` and periapsis; the belt's rms eccentricity is the quadrature sum for anything that reads it. A lone star draws today's numbers (tested exactly). Epoch with `--affects fill`; the member rows change, the belt's mass and extent do not.

## Done when

A belt in a system with a close companion shows a higher mean eccentricity than the same belt without it, by the ratio the census predicts; a belt beside a heavy giant is changed less; a lone-star system, and any system whose companions are too far to matter (`alpha` below the registered floor), gives today's numbers; stand-ins are registered in GEN.725 with burn-down lines.

## Questions for TODO (plain words)

1. The companion's periapsis direction is not stored anywhere; slice 0 draws it per system from a stream of its own. Is that acceptable, or should the star orbit store it (a column and a migration)?
2. Planets are not forced by a companion star today. File that as its own item (it would reuse slice 0)?
3. Mutual inclination of star orbits does not exist; forced inclination and Kozai pumping wait for it. File that as its own item?
