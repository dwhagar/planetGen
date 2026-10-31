# Planetary rings (GEN.732)

Status: built. The model is `src/planetgen/physics/rings.py`; the fields are stored from schema v160 (`planets` and `moons`, migration 0160) and v171 (`rogue_planets`, migration 0171); the painter (MAP.275) reads them and the colour of a ring is GEN.670's research, so neither is here.

## What a ring is here

One summed ring system per body: a kind, an origin, an inner and an outer radius (km from the planet's centre), a mass (kg), a mean optical depth over the area it covers and a water-ice share. Saturn's seven named rings and its gaps are one `broad` system from 74 500 to 140 220 km; the gaps are the painter's, from the moons' resonances. A moon draws like a planet of its own kind (GEN.836: `rings.apply` on each moon, cold when its parent is beyond the snow line, no moons of its own, so a ringed moon is rare, one in a hundred), a rogue planet has none (GEN.1057: its moons and comets are not followed), and a belt has none.

## The physics that fixes the numbers

1. **The outer edge.** Material inside the Roche limit cannot gather into a moon. For a planet of radius R and mean density rho_p and ring material of density rho_m the fluid limit is 2.456 R (rho_p / rho_m)^(1/3) (Roche 1850; Murray and Dermott 1999, recalled). Saturn's A ring ends at 136 775 km against 129 900 km for ice (ratio 1.05), Jupiter's main ring at 129 000 km against 139 000 km for rock (0.93), Neptune's Adams ring at 62 933 km against 71 300 km (0.88) and Uranus's epsilon ring at 51 149 km against 67 500 km (0.76). The outer edge is drawn as 0.75 to 1.1 of the limit, which holds all four, and held inside 0.95 of the innermost moon's distance.
2. **The inner edge.** 1.1 to 1.7 planet radii (Saturn's D ring 1.11, Jupiter's halo 1.29, Uranus's 6 ring 1.64, Neptune's Galle ring 1.7), held below the outer edge. A planet whose room between the two closes (a moon just outside the planet) has no ring.
3. **The optical depth is an output.** For particles of radius r and density rho_m, a ring of mass M over the annulus area A with a covered share f has the mean optical depth 3 M / (4 rho_m r A f). With r = 1 m, Saturn's rings (1.5e19 kg) come out at 0.27 and Uranus's (1.5e16 kg over 0.3 percent of the annulus) at 0.55, both inside what the Voyager and Cassini occultations measured.
4. **The ice.** A broad ring beyond the star's snow line keeps its ice (Saturn's are 90 to 99 percent water ice); a ring inside the snow line, dust, and the dark narrow rings have none or little. The ring's density is the volume-additive mixture of ice (917 kg/m3) and rock (2500).

## What is drawn

| Draw | Rule | Registry |
|---|---|---|
| Whether, from where, and which kind | from history (GEN.1057): a moon that fell inside the Roche limit or a comet that passed inside it, none inside 0.1 AU; the kind from the ring's mass (`dusty` under 1e12 kg, `narrow` under 3e17, else `broad`) | BC20, BC22 |
| Origin | `moon_erosion` (dust), `disrupted_moon` (narrow), `disrupted_moon` or `captured_body` (broad) | BC20 |
| Ice share | per kind, drawn always, used beyond the snow line | BC23 |
| Mass | log-uniform per kind, at most 1e-4 of the planet's mass | BC22 |
| Covered share | narrow rings 1e-3 to 1e-2 of the annulus | BC22 |
| Outer and inner edge | as above | BC21 |

Every number is a labelled stand-in in the registry (`gap-filling-plan.md` rows BC20 to BC23, listed in `stand-in-values.md`) with its tier, its tweak level and what replaces it.

## When it runs

`rings.apply(planet)` runs in the system's ring pass, after the moons and the snow line are final (`generation/system.py`), and again for a planet the validation or an admin edit gives a new class (`generation/validation.py`); `planets.generate_planet_properties` clears the old ring when a body is regenerated as another class. The pass draws from the system's stream after the life pass's predecessors, so it moved every later draw: a new generator epoch.

## Left out

Ring colour and albedo (GEN.670), the painter and the gaps (MAP.275), ring shadows and their effect on climate, the lifetime of a ring and the age at which it formed, ring moons (a moon born from a ring, `moon_composition.RING_MAX_ROCHE_RATIO`, still reads `has_ring` False: wiring it is a follow-up for GEN.466), rings of small bodies built as a module but not wired (the next section).

## Rings of small bodies at a spin-orbit resonance (GEN.1056, lock-free step)

`physics/small_body_rings.py` (registry RR1 to RR4, section 6.71) is the R4 row of weird-worlds.md section 9: a small body (under 5e22 kg) whose rotation is slow enough or whose shape pulls hard enough holds a narrow ring where the ring goes round three times for each turn of the body, past the Roche limit of its ice and inside half its Hill radius. The orbit follows from Kepler's third law (`resonance_radius_km`), 4 200 km for Quaoar (observed 4 057), 2 290 for Haumea (observed 2 287) and 426 for Chariklo (observed 391 and 405, whose mass is uncertain); the chance is a stand-in rate (0.10, no survey behind it). `ring_survives` is the R9 gate: no ring inside 0.1 AU of a star (Schlichting and Chang 2011 via the design note, 2ND-HAND). Both are pure functions that nothing calls yet: GEN.1057 wires them with the conditional rings, and the origin there is `disrupted_moon` (the cause named for Quaoar's ring), so no schema change is needed.

Ring colour and albedo (GEN.670), the painter and the gaps (MAP.275), ring shadows and their effect on climate, the lifetime of a ring and the age at which it formed, ring moons (a moon born from a ring, `moon_composition.RING_MAX_ROCHE_RATIO`, still reads `has_ring` False: wiring it is a follow-up for GEN.466), rings of small bodies (Chariklo, Haumea, Quaoar).

## Ring origin from history (GEN.1054)

`physics/ring_origin.py` is the history rule that replaced the chances of `rings.draw_ring` (GEN.1057 wired it). A ring exists when something came apart inside the Roche limit and has not spread away: a moon that fell (`moon_orbits.Fate`, `time_s`; the share `RING_RETAINED_FRACTION` of its mass stays in orbit) or a comet that passed inside the limit (a Poisson rate scaled from Saturn by the Roche sphere's cross-section with gravitational focusing, `COMET_CAPTURE_FRACTION` of it stays). Each ring spreads by 1/e every `RING_LIFETIME_YR`, none stays inside 0.1 AU of the star, the sum is held under a ten-thousandth of the planet, and the origin is that of the event holding most of it. A giant with nothing in its past has no ring. The tests are Saturn (a Mimas-sized moon 150 Myr ago leaves 1.5e19 kg to a factor of two) and Phobos. The stand-ins of the module are RO1 to RO8 (plan 6.70).

## Built: rings from history (GEN.1057)

`rings.apply(planet, host)` no longer draws whether a planet has a ring. `moon_history` keeps the moons the tide pulled inside the Roche limit as `planet.fallen_moons` (mass and the year they fell; a moon dragged below the innermost orbit counts at the planet's age), `ring_origin.draw_comet_passages` draws the comets that passed inside the limit over the star's age, and `ring_origin.ring_origin` sums their rings, each spread by 1/e every half billion years, capped at a ten-thousandth of the planet and gone inside 0.1 AU. The ring's origin is the event that holds most of it (`disrupted_moon` or `captured_body`; `moon_erosion` is no longer drawn), its age that event's (`ring_age_yr`, schema v187, the years since), its kind from its mass (`rings.kind_of_mass`) and its edges and ice share from the old draws of `rings.build_ring`. A dwarf planet that orbits a star and has no such ring may hold the narrow ring of its spin-orbit resonance (`small_body_rings.resonance_ring`, R4); it has no age. A moon has the ring its history made (a moon has no moon to lose and a tiny Roche sphere to catch a comet, so almost never one); a rogue planet has none, because nothing follows its moons or its comets. Ice share and darkness by age and pollution (R6, GEN.1055) and the ring-plane follow-up of the spin (R3, GEN.1058) are not built.

**Brightness and shadow (GEN.1055).** `physics/ring_darkening.py` is the pure module for how bright a ring is and how much of the sunlight it shades: meteoroid dust darkens the ice at a rate set by the ring's mass per area, and the rings shade the latitudes behind them by `exp(-tau/sin(declination))` over the lit half of the day. `rings.apply` calls it (`rings.darken`) to store the ring's ice share as it is now; the climate does not read the shadow yet.
