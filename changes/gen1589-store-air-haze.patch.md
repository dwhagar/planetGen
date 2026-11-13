**Headline:** A world whose air life has changed now keeps the temperature its haze gives it, and its carbon dioxide is solved against the carbon sink to the end.

### Changed
- **`store_air` leaves out the haze (GEN.1589).** The air a life stage stores (the oxygen rise, the carbon balance, the trace gases) was warmed as clear air, while the air solver had put a share of the sunlight into a haze high up (`aerosols.haze_of`). A body now keeps its haze draws, its star's temperature and ultraviolet and the heat from inside (`_air_haze`, runtime only), and `biogenic_oxygen.store_air` works the haze and the temperature against each other (`HAZE_PASSES`, damped) as the air solver does. Of the 187 bodies in 3 x 40 seeds that a life stage stores, 79 had a stored temperature that a fresh solve of their stored air gave differently by more than 1 K (32 by more than 20 K); with the haze 61 of 191 do and 21 by more than 20 K. The rest are the 2000 K clamp and the swings of the greenhouse table for a trace gas (GEN.1588).
- **`rebalance_co2` stopped before it settled (GEN.1589).** Its damped passes stopped after 8 with the answer still falling (562 K to 432 K in one air). The CO2 the air asks for less the CO2 the sink holds is a rising difference of the CO2, so the root is found by bracketing (doubling up to the retention cap) and halving to `REBALANCE_TOLERANCE`; a run that does not settle is logged.

generation-output: changed
generation-affects: fill
