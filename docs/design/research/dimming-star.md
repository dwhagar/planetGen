# Design note: the dimming star (GEN.1437, part 2)

Written by Lane 21, 2026-11-11, from section 5.3 and item A7 of `/mnt/project-files/reports/flavor-text-and-sci-fi-anomalies-2026-11-08.md`. Design only: nothing here is built. The fictional organism of that item is out of scope; the physical version is a star whose light dips irregularly because dust crosses our line of sight to it. Numbers not tied to a reference are proposals and stand-ins to be registered (GEN.725) when built.

## 1. What the real cases show

- **KIC 8462852 (Tabby's star).** Aperiodic dips of up to about 20 percent lasting days to weeks, in the Kepler data of 2009 to 2013 ([R5220], Boyajian and others 2016, abstract recalled; the paper's tables not re-read), and more dips after Kepler with a colour that points to dust ([R5221], Boyajian and others 2018, abstract recalled). It was one star in about 150,000 Kepler targets, so the rate is at most about 1e-5 per star; the cause is not settled (dust from a disrupted body or a comet swarm are the physical candidates).
- **Disintegrating planets.** Kepler-1520 b ([R3666], Rappaport and others 2012), K2-22 b ([R3903], Sanchis-Ojeda and others 2015) and KOI-2700 b: a rocky planet losing dust in a tail, with a transit depth that changes from transit to transit (about 0 to 1.3 percent, recalled). K2-22 b's grains were measured at a large size and a high loss rate by [R3954]. The generator has the cause as the oddity W14 (`physics/disintegrating_planet.py`: the wind, the mass loss and the lifetime, from [R3487]); the module's own text says "the tail and the transit depth are not here".
- **Young stars behind an inclined disk (UX Ori type).** A pre-main-sequence star of 2 to 8 solar masses (a Herbig Ae/Be star) or a T Tauri star dims when the clumps of its inner disk cross the line of sight. The generator already names these stars (`physics/variability.py`: `T_TAURI`, `HERBIG`) and shows a Variability row; it has no dimming for them.

## 2. What the generator already has, and does not

Has: the belts and their dust (`physics/debris_dust.py`: fractional luminosity of a belt, GEN.1071); the disintegrating planet's mass-loss rate (above); fallen moons and their rings (`moon_history`, `ring_origin`); the variable-star kinds (above); and a star's disk while it is young. Does not have: any observer. The generator draws a system without a line of sight, so a dip is not something the system "has" until it is seen from somewhere. That fixes the design: **store the cause and the largest dip a sight along the dust's plane would show**, not a light curve.

## 3. The design: a dimming state with three causes

A star has a **dimming state** when a source of dust big enough to matter lies in or near the plane that would cross it. The state is derived from stored values; nothing is rolled except where a cause already is a draw.

1. **Dust stream of a disintegrating planet** (W14). The tail's cross-section, from the mass-loss rate and a grain size (the grain size is a stand-in; [R3954] measured it for K2-22 b). The dip: the tail's covering fraction of the stellar disk, which is the optical depth times the tail's width over the star's diameter. Periodic with the planet's period; depth varies because the loss is not steady (the planet's `DUTY_CYCLE`, `disintegrating_planet`). Largest dip near 1 percent for a planet of the Kepler-1520 b type.
2. **Occulting disk of a young star.** For a star in the `T_TAURI` or `HERBIG` kind, the inner rim of its disk is at the dust sublimation radius (`planet_architecture.sublimation_au`, already used for the orbit floor), and the disk's scale height over the radius sets how many degrees from the plane a sight is still blocked (a stand-in: aspect ratio 0.05 to 0.1). The dip is the clump's optical depth; the state says "dims irregularly when seen within N degrees of the disk plane".
3. **Dust from a recent breakup** (the Tabby case). A belt or a ring whose dust is far above what the belt's grinding gives (a recent collision or a disrupted moon or comet swarm) crossing the star. Its trigger is a draw already in the chain: the stored belt dust (`debris_dust`) against a threshold, or a ring from a disrupted moon (R1 of the weird-worlds table) inside the Roche limit with the star at the dust's distance, and a rate no higher than the 1e-5 per star above. The dip: the dust's covering fraction. This is a case where the cause is not known in the real star, so the row is a stand-in of the lowest tier (from needed output) and says so in plain words.

For all three: the stored fields are `dimming_kind` (one of the above or empty), `dimming_max_depth` (the largest fraction of the star's light the dust can block when seen edge-on) and `dimming_plane_deg` (how close to the dust's plane a sight must be, from the aspect ratio or the tail's width over its distance). A page reads them as "dims by up to X percent, irregularly, if seen within Y degrees of the plane of the dust". The sensor reading in the fiction (an infrared-bright arc between star and planet) is the dust's infrared excess, which `debris_dust` can already give as the fractional luminosity.

## 4. How we will check it

Limits: `dimming_max_depth` under 1 for any star, and under 0.02 for the planet-tail cause (the observed depths); the plane angle under 30 degrees for the disk cause; a dimming state only for a star with a source of dust; the tail's mass loss times the age under the planet's mass (the W14 limit); the share of stars with a state under the rate stated. Calibrate the three covering fractions against the depths above (Kepler-1520 b about 1 percent; Tabby's up to about 20 percent as an upper bound for the breakup cause).

## 5. Chain items this needs (not built here; TODO numbers them)

1. **Tail model for W14**: the dust tail's optical depth, width and the transit depth, from the loss rate and a grain size, with the stand-in rows. Lock-free (a pure function and tests on Kepler-1520 b and K2-22 b).
2. **Young-star disk occultation**: the aspect-ratio rule and the plane angle, with tests on a Herbig and a T Tauri star. Lock-free.
3. **Breakup dust**: the trigger from the belt dust and from a ring of fallen moons, the rate cap, the stand-in row at the lowest tier. Lock-free pure function.
4. **Stored fields** `dimming_kind`, `dimming_max_depth`, `dimming_plane_deg`, with columns, sheet rows, unit registry entries, and a Variability row line (an epoch and a migration).
5. **The oddity row** (a gate on the stored `dimming_kind`, no extra rate), the gated flavor lines, and the plausibility limits of section 4.
