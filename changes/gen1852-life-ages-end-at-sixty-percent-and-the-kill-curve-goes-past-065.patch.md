**Headline:** A life age now ends at an extinction that kills 60 percent of species (the Earth's Big Five rate), and the kill curve goes on past 0.65 to the end-Permian's 0.90.

### Changed

- GEN.1852: `life_ages.AGE_KILL_THRESHOLD` is 0.60, the Boss's pick of 11:18 PM EST 11-26 (it was 0.30). On 227 generated systems (1,384 bodies with life, 4,750 body-Gyr) the 0.30 threshold gave 23.7 events a Gyr against about 9 for the Earth's Big Five; 0.60 gives 9.7. The impact that reaches it is the one the Earth has once in 72 Myr (about 10 km, was 4.1 km), so a body now draws fewer, larger impacts for its life ages and its extinctions (`air_events.events_of` is called with the larger smallest impact), and the flood basalts (0.46 on the curve) no longer end an age. The extinctions list and the "extinction floor" on the page follow the constant.
- GEN.1875: the kill curve (`air_events.KILL_CURVE`) ended at 0.65 for 100 Myr and held it for every rarer event, a cap with no source behind it (the table's last point held flat; the two Raup abstracts reachable give nothing past 100 Myr). It now has a fourth point, the end-Permian: 0.90 at 540 Myr (once in the Phanerozoic, inside the 0.81 to 0.96 the estimates give), straight in the log from 0.65 at 100 Myr, held at 0.90 beyond. Registered as EV6b (tier T3, tweak L1). Impacts over about 12 km now kill between 0.65 and 0.90 of species where they all killed 0.65; nothing under it moves.

generation-output: changed
generation-affects: fill
