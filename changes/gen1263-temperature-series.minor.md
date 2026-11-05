# GEN.1263: the temperature of every place at every time

### Changed

**Headline:** A new calculation gives the temperature of every place on a world, by latitude, longitude and time of day, through its whole year (no stored value changes yet).

- `seasonal_swing.temperature_series` extends the seasonal energy balance with the longitude and the time of day on the 12 by 12 cells of `surface_cycle`. One calculation serves every body: a tidally locked world, a world that turns backwards, a slow spinner near resonance and Earth.
- The heat capacity follows the cycle: the air column, the sea to the depth the cycle mixes (a few metres for a day, 50 m for a year) on the ocean share, and the ground's skin on the rest. The air carries heat toward the world's mean (`transport_share`, the locked share for a body that turns once an orbit).
- Ice is laid by place from the stored ice share with the GEN.972 law, so cold cells are brighter.
- The series repeats exactly (`surface_cycle.periodic_span`), the outgoing emission is scaled until its area-and-time mean is the stored mean, and a Fourier filter with Anderson mixing solves the fourth-power emission.
- Eight stand-ins are registered in GEN.725 (rows SS4 to SS11 of section 6.91). Nothing calls the series yet (GEN.1264 reads shares and extremes off it; the score, columns and migration are GEN.1265 and GEN.1266), so no generated output changes.
