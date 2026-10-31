**Headline:** The Search page now has range boxes for gravity, temperature, pressure, radiation dose and the other physical values, each with its unit beside it, and shows the ranges as chips.

### Added
- **Search by physical values.** A new fold on the Search page, with a row for each value of planets, moons and stars: surface gravity, surface temperature, surface pressure, surface radiation dose, orbit distance, day and year length, rotation period, mass, magnetic moment, ocean depth, number of moons, and for stars luminosity, age and temperature. Each row has a least and a most box and a unit menu (for example bar, kPa or atm for pressure). Typing 0.8 to 1.2 g and a dose of at most 2 mSv a year finds the planets that have both.
- **Gases.** Rows for the common gases (nitrogen, oxygen, carbon dioxide, water, methane, argon, hydrogen, helium, ammonia, sulfur dioxide) give the partial pressure at the surface; a row called "another gas" reaches any other gas the generator knows.
- **Chips and a Values column.** Each range is a chip with its unit and a cross that takes only that range away. The planet, moon and star tables have a new Values column that shows the values the search asked for, with their units.
- **Order by a value.** A menu for each kind of object orders its results by any of the values, smallest or largest first; a table header's own sort still wins.

### Changed
- **The address stays short.** A unit that is the value's own, or that has no range beside it, is not put in the address.
