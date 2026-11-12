**Headline:** a planet or moon is now classed after its latitude bands exist, so the classes that read the hottest and coldest region, the night side and the wet bulb (M, O, H, P and R) can match when a body is built, and a moon is classed after its water is remade at the final age.

### Changed

- `planets.apply_class` is split into `prepare_class` (the body type and the roundness ratio, which the spin chain reads) and `classify` (the letter, its match, the runner-up and the digit). The build (`Planet.__init__`, `finish_climate`), the rebuilds (`remake_end_stages`, the giant settle, the late-accretion remake) and the spacing pass (`validation`) now prepare the class, make the spin (whose last step, `habitability_world.generate`, makes the bands), and classify after. `apply_class` stays for callers that already hold the bands.
- A moon's class is applied in `_finish_moon`, after `refresh_at_final_age` remakes its hydrosphere, aquifer and groundwater, and not in `_finish_planet` before it.
- A request for a class (`FORCE_CLASS`, a slot, the habitable-world request) now sees the class of the finished conditions, so the draws it compares are the draws that can match. Whether M can match at all is the Class M audit's (GEN.1581, GEN.1582, GEN.1585, GEN.1586); this change removes the ordering that made it impossible.
- Changes the class letter of any body whose class read the bands, so it is a lock item: the golden and epoch are the train's.

generation-output: changed
generation-affects: fill, population
