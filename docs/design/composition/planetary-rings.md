# Planetary rings (GEN.732)

Status: built. The model is `src/planetgen/physics/rings.py`; the fields are stored from schema v160 (`planets` and `moons`, migration 0160) and v171 (`rogue_planets`, migration 0171); the painter (MAP.275) reads them and the colour of a ring is GEN.670's research, so neither is here.

## What a ring is here

One summed ring system per body: a kind, an origin, an inner and an outer radius (km from the planet's centre), a mass (kg), a mean optical depth over the area it covers and a water-ice share. Saturn's seven named rings and its gaps are one `broad` system from 74 500 to 140 220 km; the gaps are the painter's, from the moons' resonances. A moon draws like a planet of its own kind (GEN.836: `rings.apply` on each moon, cold when its parent is beyond the snow line, no moons of its own, so a ringed moon is rare, one in a hundred), a rogue planet draws with no host and no moon (`rings.apply_rogue`, cold below the snow-line temperature of 170 K, a brown dwarf has none), and a belt has none.

## The physics that fixes the numbers

1. **The outer edge.** Material inside the Roche limit cannot gather into a moon. For a planet of radius R and mean density rho_p and ring material of density rho_m the fluid limit is 2.456 R (rho_p / rho_m)^(1/3) (Roche 1850; Murray and Dermott 1999, recalled). Saturn's A ring ends at 136 775 km against 129 900 km for ice (ratio 1.05), Jupiter's main ring at 129 000 km against 139 000 km for rock (0.93), Neptune's Adams ring at 62 933 km against 71 300 km (0.88) and Uranus's epsilon ring at 51 149 km against 67 500 km (0.76). The outer edge is drawn as 0.75 to 1.1 of the limit, which holds all four, and held inside 0.95 of the innermost moon's distance.
2. **The inner edge.** 1.1 to 1.7 planet radii (Saturn's D ring 1.11, Jupiter's halo 1.29, Uranus's 6 ring 1.64, Neptune's Galle ring 1.7), held below the outer edge. A planet whose room between the two closes (a moon just outside the planet) has no ring.
3. **The optical depth is an output.** For particles of radius r and density rho_m, a ring of mass M over the annulus area A with a covered share f has the mean optical depth 3 M / (4 rho_m r A f). With r = 1 m, Saturn's rings (1.5e19 kg) come out at 0.27 and Uranus's (1.5e16 kg over 0.3 percent of the annulus) at 0.55, both inside what the Voyager and Cassini occultations measured.
4. **The ice.** A broad ring beyond the star's snow line keeps its ice (Saturn's are 90 to 99 percent water ice); a ring inside the snow line, dust, and the dark narrow rings have none or little. The ring's density is the volume-additive mixture of ice (917 kg/m3) and rock (2500).

## What is drawn

| Draw | Rule | Registry |
|---|---|---|
| Whether and which kind | a giant always has a ring, `dusty` : `narrow` : `broad` = 1 : 2 : 1 (Jupiter; Uranus, Neptune; Saturn); a terrestrial or dwarf planet has one with probability 0.01 | BC20 |
| Origin | `moon_erosion` (dust), `disrupted_moon` (narrow), `disrupted_moon` or `captured_body` (broad) | BC20 |
| Ice share | per kind, drawn always, used beyond the snow line | BC23 |
| Mass | log-uniform per kind, at most 1e-4 of the planet's mass | BC22 |
| Covered share | narrow rings 1e-3 to 1e-2 of the annulus | BC22 |
| Outer and inner edge | as above | BC21 |

Every number is a labelled stand-in in the registry (`gap-filling-plan.md` rows BC20 to BC23, listed in `stand-in-values.md`) with its tier, its tweak level and what replaces it.

## When it runs

`rings.apply(planet)` runs in the system's ring pass, after the moons and the snow line are final (`generation/system.py`), and again for a planet the validation or an admin edit gives a new class (`generation/validation.py`); `planets.generate_planet_properties` clears the old ring when a body is regenerated as another class. The pass draws from the system's stream after the life pass's predecessors, so it moved every later draw: a new generator epoch.

## Left out

Ring colour and albedo (GEN.670), the painter and the gaps (MAP.275), ring shadows and their effect on climate, the lifetime of a ring and the age at which it formed, ring moons (a moon born from a ring, `moon_composition.RING_MAX_ROCHE_RATIO`, still reads `has_ring` False: wiring it is a follow-up for GEN.466), rings of small bodies (Chariklo, Haumea, Quaoar).

## Ring origin from history (GEN.1054)

`physics/ring_origin.py` is the history rule that will replace the chances of `rings.draw_ring` (GEN.1057 wires it). A ring exists when something came apart inside the Roche limit and has not spread away: a moon that fell (`moon_orbits.Fate`, `time_s`; the share `RING_RETAINED_FRACTION` of its mass stays in orbit) or a comet that passed inside the limit (a Poisson rate scaled from Saturn by the Roche sphere's cross-section with gravitational focusing, `COMET_CAPTURE_FRACTION` of it stays). Each ring spreads by 1/e every `RING_LIFETIME_YR`, none stays inside 0.1 AU of the star, the sum is held under a ten-thousandth of the planet, and the origin is that of the event holding most of it. A giant with nothing in its past has no ring. The tests are Saturn (a Mimas-sized moon 150 Myr ago leaves 1.5e19 kg to a factor of two) and Phobos. BC20 to BC23 stay as the draw until then and are T4, L3; the stand-ins of the module are RO1 to RO8 (plan 6.69).
