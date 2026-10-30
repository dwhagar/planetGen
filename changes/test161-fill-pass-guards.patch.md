**Headline:** New tests hold each kind's system chance to its table, and keep planets and moons out of a scatter row before its sector fills.

### Changed

- TEST.161: tests/test_fill_pass_guards.py checks that the share of every kind's objects that roll a system follows the host tables within four Poisson sigma (a kind the table gives nothing never rolls one), that a rogue planet keeps its moon at the table's chance, that a scatter row has no column for a planet, moon or belt, that a scatter stores its rows and writes no planet, moon, belt or rogue planet before its sector fills, that a filled sector with a black hole holds no planet of its own, and that a body past the orbit ceiling goes outermost first. The stream independence and same-seed same-system checks already in test_gen_fill_stages.py cover the other two items. Tests only: no generated value changes.
