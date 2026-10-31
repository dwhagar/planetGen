**Headline:** Moons now drift, fall into rings or are lost under the tide of their planet over the age of the system, and the tides they raise heat the planet.

### Changed
- **The moons' history (GEN.959, stage 3 of the tides chain).** Each moon's distance from `generate_moons` is where it is born; the tide it and the star raise on the planet then moves it for the age of the star (`moon_history`). A planet that turns faster than its moons go round pushes them outward and slows (the Earth and the Moon); one the star has locked, or a slow one, pulls its moons in. A moon that crosses its Roche limit is removed and its mass becomes a ring of disrupted moon on the planet; one pushed past the stable part of the Hill sphere is lost. About two thirds of the moons of planets under 0.3 Earth masses and a quarter of the moons of larger planets are gone over 40 default systems.
- **The planet's day follows the exchange.** A planet with moons ends at the spin the tides leave, or locked to its star when it ends within a thousandth of the star's orbital rate; the thermal tide of its air (GEN.957) is in the same sum.
- **The moons heat the planet.** A rocky planet's tidal heat now includes the power the tides of its moons dissipate in it, at the mantle's `k2/Q` for each tide's frequency (an Earth under the Moon: 2 TW), in the heat flow, the melt and the heat pipe.
- **Stand-ins.** TM5 to TM9 (plan section 6.25).

generation-output: changed
