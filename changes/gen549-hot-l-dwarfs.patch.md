**Headline:** Brown dwarfs near the hydrogen-burning limit stay hot, so L dwarfs hotter than 1,500 K are as common as the 20 pc census counts them.

### Changed
- A brown dwarf's temperature and radius at an age come from the cooling track of Baraffe et al. (2015), continued past a track's last age at Burrows et al.'s old-age exponent, in place of a single power law that was too cool for the heaviest brown dwarfs. Its luminosity is the Stefan-Boltzmann value of the two.
- A free-floating brown dwarf draws its age first, then takes its radius, temperature and luminosity from the track at that age, so the three agree.
- `tests/test_brown_dwarf.py`: the expected failure for L dwarfs hotter than 1,500 K is gone; they now reach the census's share within a factor of 2. The design note `brown-dwarf-state.md` names the track and its limits.

generation-output: changed

TODO items resolved: GEN.549 / Added: none
