**Headline:** A giant that swept an outer disk now holds a trojan swarm at its L4 and L5, as two regions that are stored, named for the giant and shown on their own page.

### Added
- GEN.1150 slice 1 (lock item, migration 0269, number settled at merge): the mass the scatter step gives a giant's 1:1 zone (`swept_hand_on.trojan_part`: the bound share of what its outer disk lost, times the swept split's `trojans` row, times what the other giants and companion stars leave it, GEN.1455) becomes a leading (L4) and a trailing (L5) `TrojanRegion` (`generation/trojan_region.py`), drawn by the existing `trojan_swarm.draw_regions` (the split, each region's libration, mean eccentricity and inclination) on streams of the swarm's own. A region is a belt-like row of `asteroid_belts` whose parent is the giant, with the new columns `lagrange_point`, `libration_amplitude_deg`, `libration_phase_deg`, `mean_eccentricity` and `mean_inclination_deg`. It is made after the spacing pass, the irregular moons and the comets, so no pass sees a body that shares a giant's orbit, and it is kept only above the belt minimum mass. Its name is `<giant> L4 trojan region`, worked out and never stored, and it takes no number among the star's belts.
- No rate is set. In 160 G2V systems 19 held a giant and 2 had an outer disk it had swept (the swarm needs one, and a giant of a high Safronov number throws nearly all of it out): about 1 percent of G2V systems hold a swarm worth keeping, the heavier of the two 1.1e22 kg. Registered as GEN.725 row TJ12.
- Not yet: the libration overlay of the orbit updater and the map icon (slice 2, MAP.354), members of a region, and whole planets at L4 and L5 (GEN.1151).

generation-output: changed
generation-affects: fill
