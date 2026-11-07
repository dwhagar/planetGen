**Headline:** A living world's oxygen is now held to a fire limit of oxygen pressure (35.5 kPa, 35 percent of one Earth atmosphere) and not to 35 percent of its air, so a thin air can carry a larger share of oxygen and the pile-up of worlds at exactly 35 percent is gone (GEN.1320).

### Changed

- `biogenic_oxygen._fire_limit_mol` solves the oxygen the other gases leave room for under a partial-pressure cap, with the share held to 95 percent; `mixture_oxygen.BIOTIC_O2_MAX` is replaced by `BIOTIC_O2_PRESSURE_MAX_PA` (stand-in OX1, section 6.195).
- Thick airs are unchanged in practice (35.5 kPa of oxygen in 1.2 MPa is still a small share); the class windows that read the share of oxygen are judged by pressure in GEN.1321.
