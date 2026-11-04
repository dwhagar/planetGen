# Day and night time shares and the hours of daylight (GEN.767)

Boss, 7:54 PM and 8:22 PM EDT 2026-10-24: every object that has a day or night gets its day and night time shares, a moon with the
extra night its planet's shadow adds, and the planet says how many hours of daylight it has. The module is
`src/planetgen/physics/daylight.py`; it draws nothing and has no stand-in value, it reads the stored rotation period, orbital
period, spin axis, orbit planes and the star's radius.

## 1. What is stored

| Column | Meaning |
|---|---|
| `day_share_equator` | share of the time the star is up at the equator, averaged over a circular year (`NULL` when the world turns once an orbit) |
| `day_share_mean` | the same averaged over the whole surface |
| `parent_shadow_share` | a moon only: share of its time in its planet's shadow |
| `solar_day_hours` | hours from one noon to the next (`NULL` when the world turns once an orbit) |
| `daylight_hours` | the solar day times `day_share_equator` |
| `polar_circle_deg` | latitude beyond which the star stays up or down for at least a day |

The night share is one minus the day share.

## 2. The rules

* The star is up while its upper limb is above the horizon. Its angular radius is `alpha = asin(R / d)`; at a latitude `phi` and a
  declination `delta` the half day arc has `cos H = (-sin alpha - sin phi sin delta) / (cos phi cos delta)` and the share is `H / pi`.
* The year is circular (the model has no eccentricity), so `sin delta = sin e sin lambda`, `e` the angle between the spin axis and
  the normal of the orbit round the star (a moon's: its spin axis against its planet's orbit normal). A tilt past 90 degrees is a
  retrograde spin and gives the same seasons.
* The surface mean is exact: a sphere is lit over `(1 + sin alpha) / 2` of its area at every instant, whatever the tilt.
* The solar day is `1 / |1/P_rot -/+ 1/P_orb|`. A world with `P_rot = P_orb` has no day; its day side is lit all the time, the
  terminator half of it, the night side never.
* The polar circle is `90 - e - alpha` degrees.
* A moon is in its planet's shadow while the centre of the star is behind the planet's disc as the moon sees it
  (`theta = asin(R_planet / a)`). With `beta` the angle between the star and the moon's orbit plane that is a share
  `arccos(cos theta / cos beta) / pi` of an orbit when `beta < theta`, and the year turns `beta` through `sin beta = sin(Theta)
  |sin lambda|`, `Theta` being the angle between the moon's orbit normal and its planet's. Eclipses come in seasons: an inclined
  orbit has none outside them, but the star crosses its plane twice a year, so no circular orbit is free of them.
  Both day shares of the moon are multiplied by one minus that share.

## 3. Checks (hand work, in `tests/test_daylight.py`)

Earth: 50.15 percent at the equator, 12.04 h of daylight in a 24 h solar day, polar circle at 66.3 degrees. Venus's solar day is
116.75 Earth days. Io in Jupiter's shadow 5.3 percent of its orbit (2.2 h of 42.5 h); the Moon 0.05 percent of the time (about four
hours a year, in eclipse seasons of about 43 days).

## 4. Stand-ins and limits

None are values. Limits: circular orbits; the star is a disc of its radius at the parent's distance (a circumbinary host uses its
primary's radius when it has one, otherwise a point source); no refraction; a moon locked to its planet loses the extra night on its
near side only, which the mean does not show.
