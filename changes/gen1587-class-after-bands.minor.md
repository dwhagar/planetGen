**Headline:** A planet or moon is now classed after its latitude bands exist and after its water is made again, so requests for a class that reads the bands (M, O, H, P, R) can be met and a moon's class follows its final water.

### Fixed
- **GEN.1587: the class was judged before the bands.** `apply_class` ran before the spin (`generate_spin` ends by making the bands through `habitability_world.generate`), so M, O, H, P and R, which read the hottest and coldest region, the night side and the wet bulb, scored 0 at build: a request for one of them (`FORCE_CLASS`, a slot) threw every draw away and was left unmet. `physics.planets.apply_class` is now `prepare_class` (the type and the roundness ratio) and `classify`; a planet is prepared before its spin and classified after it in `Planet`, in `finish_climate`, in `reconcile_moved_planet`, in `remake_end_stages` and when the air is made again. No second class pass: the call that sat before the spin moved after it.
- **A moon's class was judged on the water it was drawn with.** `_finish_planet` classed each moon before `_finish_moon` made its hydrosphere, aquifer and groundwater again at the final age. The moon is now classed in `_finish_moon`, after `refresh_at_final_age`.
- Changes the class of bodies (and so which draws a request keeps) wherever a band-reading class applies, so it is a lock item: the golden and the epoch are the train's.

generation-output: changed
generation-affects: fill, population
