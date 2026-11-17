**Headline:** A planet that loses all its moons after they were made no longer keeps the wobble they gave it.

### Fixed
- A planet's reflex offset (its own wobble from its moons) was made when the moons were and left alone when the planet had none, so a planet whose moons all went afterwards (a Kozai-Lidov loss in `_kozai_moons`, moons that no longer fit after its class changed, `drop_unstable_moons`, a bar on habitable worlds) still carried the old offset, for example 3.0e-5 AU (an Earth-Moon pair is about 3e-5 AU) on a moonless planet. `StarSystem._refresh_wobble` and the two refreshes in `validation` now make it again from the moons that are left, and zero it with none. A test of `_refresh_wobble` holds it; `test_orbital_motion` generates each of its systems under a stated seed (`orbital-motion:N`) so a failure replays.

generation-affects: none

generation-output: unchanged
