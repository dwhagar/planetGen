**Headline:** A rogue planet thrown out of a young system may now leave bound to a companion, at the share seen in young clusters, and its page names the distance and the period of the pair.

### Added

- **Rogue planets that stay a pair (GEN.332, slice D).** Free-floating planetary-mass objects are seen in pairs: about 9 percent of the systems of free planets in the Trapezium sample are bound, 28 to 384 AU apart (`physics/rogue_pairs.py`, [R3467], read second hand). The pair does not last, so the share falls by a factor of 50 in every 10 Myr after the ejection (Portegies Zwart and Hochart 2024, second hand): the rogues of a young donor may have a companion, a rogue thrown out a hundred million years ago or more never does. `rogue_yield.bind_pairs` decides it from the donor's own stream (the donor's rogues stay a function of the donor alone): a body not yet taken is a pair with the chance of its age, its companion is the next body the donor ejected from the same star by the same cause, which leaves with the same moment and kick and lands with it. The separation is log-uniform in the observed range and held inside what the galaxy's tide leaves a pair, the period is Kepler's third law on a circular orbit (no source gives an eccentricity for wide young pairs). The rogue's note says it left bound to a companion, and the rogue planet page and sheet give the distance to it and the period.
- **Storage.** Migration 0260 (schema v260) adds `rogue_planets.pair_event_index` (the companion's event index, the link both rows carry), `pair_separation_au` and `pair_period_yr`, all NULL for a rogue that left alone and for every row made before v260.
- Stand-ins RPP1 to RPP6 in section 6.290 of the gap-filling plan and the list of every stand-in (GEN.725): the share (tier T3, tweak L1), the survival (T3, L0), the separation range (T3, L0), the tidal limit (T5, L0), the circular orbit (T7, L1) and the pairing by order (T6, L1).

generation-output: changed
generation-affects: fill
