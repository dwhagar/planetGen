**Headline:** The search can now keep and order planets, moons and stars by their stored physical values, with the bounds given in the unit you choose.

### Added
- **Range filters on physical values.** `GET /api/search` takes `planet_min_gravity=0.8&planet_max_gravity=1.2&planet_max_surface_dose=2` (and the same with `moon_` or `star_` in front) to keep the bodies whose gravity is 0.8 to 1.2 g and whose surface dose is under 2 mSv a year. The keys are gravity, surface temperature, surface pressure, surface dose, magnetic moment, ocean depth, orbit distance, day length, rotation period, year length, mass, the number of moons, the partial pressure of each atmosphere gas, and for stars luminosity, age and temperature. A bound may be given in any unit of its kind with `<entity>_<key>_unit` (pressure in kPa instead of bar).
- **Sorts on the same values.** `planet_sort=gravity` orders a panel smallest first, `-gravity` largest first; a body with no value comes last.
- **The answer says what it used.** Each row carries the values of the keys the request used, with their units, and the answer lists every key (`values`) and the ranges applied (`active_ranges`).

### Changed
- **Search parameters are checked.** A bound that is not a number, a minimum above the maximum or a unit of the wrong kind answers `400` with the reason.
