**Headline:** The temperature of every place on a body now costs about a twentieth of what it did at full resolution, so it can run for every rocky body in bulk generation and in a recompute (no stored value changes).

### Changed

- `seasonal_swing.temperature_series` takes a default resolution of 12 rows of latitude by 4 columns of longitude (12 for a body that turns fewer than 24 times an orbit, where longitude is the place in the day), 8 steps a solar day, 120 steps an orbit and at most 1500 steps, stopping when no temperature moves more than 1 K; `FINE` is the old resolution.
- Measured on 310 rocky planets and moons drawn from four star types and three zones: a mean of 66 ms, a median of 26 ms, 158 ms at the 90th percentile and 1.0 s at the most (series and read-offs together), against a mean of 1.1 s at `FINE`. The livable share is within 0.005 of the fine one on average and 0.033 at the 95th percentile. The slow bodies are airless ones with cold dry nights, which need many passes of the fourth-power emission.
- A test holds the work of a body to its budget of points and passes and the livable share of the default to the fine one's on the ten example worlds. Registry row SS11 of GEN.725 has the new values.
