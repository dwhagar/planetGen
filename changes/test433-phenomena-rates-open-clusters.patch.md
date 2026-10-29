**Headline:** The prevalence suite now also checks the rogue planet, brown dwarf, runaway star and interstellar debris rates and the number of open clusters against the surveys.

### Added
- **The phenomena census, second cut (TEST.433).** `src/tests/audit/prevalence_phenomena.py` adds seven measures to `src/tests/data/prevalence_phenomena_expectations.json`: terrestrial rogue planets per star (2 to 10) and per pc^3 (0.5 to 1.4), Jupiter-mass rogue planets per star (at most 0.25, Mróz et al. 2017), brown dwarfs per pc^3 (0.025 to 0.035, Kirkpatrick et al. 2021), runaway stars per star (1 to 2 percent, Tauris 2015), interstellar debris per pc^3 (7e14 to 1.8e15) and the open clusters of the reference galaxy (1,867 catalogued as a floor, up to the generator's cap). The rates are read where the sectors' fill reads them, so a dialled `PHENOMENON_RATE_SCALE` trips a flag. All seven hold today.

TODO items resolved: none (TEST.433 stays open for the nebula classes). Added: none.
