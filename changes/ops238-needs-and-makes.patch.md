**Headline:** every command now says what it needs and what it makes, and a test holds that to the order a galaxy run has today.

### Added

- OPS.238, first part (the declarations): `queue/needs.py` names the artifacts of a generation run (the skeleton, the sectors, the bright-star passes, the phenomena, the nebula field, the population and the map totals) and gives each command its `needs` and `makes`. `tests/test_needs_and_makes.py` checks that every artifact has a maker, that the graph has no cycle, that the fixed stage order of a galaxy run and of each scatter run meets every need a command of the run can meet, and that everything `commands.REQUIRES` holds is a need the declarations state. Nothing reads the declarations to start a run yet: `REQUIRES` goes when the planner (the second part) reads them.
