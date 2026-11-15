**Headline:** A cool giant with a carbon-rich envelope no longer stops a galaxy run with "an abundance must not be negative".

### Fixed

- The clouds of a giant's layer (`mixture_giant.remove_clouds`) took the oxygen the oxides bind off a budget whose carbon and oxygen had been re-split to a high carbon to oxygen ratio (`element_budget`), so a metal-rich envelope at a cool layer was left a negative oxygen and the equilibrium refused it (`test_db_column_round_trip`, Train 31). The oxides now bind at most the oxygen the gas has. It showed only once GEN.1717 gave such a giant its true, cooler 1 bar temperature (the table had given it thousands of kelvin, above every condensation). No value changes where the oxygen was enough.
