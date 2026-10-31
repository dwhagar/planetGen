**Headline:** The temperate-ocean table by class now counts a surface ocean and an ice-covered ocean in separate columns.

### Changed
- **`src/tests/temperate_ocean_by_class.py` splits the ocean (GEN.1144 follow-up).** The first version counted any body with an ocean fraction over zero and a hydrosphere other than ice, so the glaciated worlds of class P (an ocean under an ice shell, whose ocean fraction is the ice cover) showed as 27 of 32 wet while `scripts/class_counts.py` found none with a surface ocean. The table now has a column for a surface ocean (hydrosphere surface ocean, hycean or brine seep, the set `class_counts.py` uses) and one for an ice-covered ocean. A tool only: no generated value changes.
