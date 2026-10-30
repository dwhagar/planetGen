**Headline:** A new audit script counts the Class D bodies, small moons, dwarf-size bodies and belts the generator makes per system, for each census star type.

### Changed

- GEN.1087: `tests/audit/dwarf_census.py` (`python3 -m tests.audit.dwarf_census --count 100`) draws full systems per star type and reports means and shares of Class D bodies (planets and moons), moons of 50 to 500 km, bodies of 400 to 2,500 km, belts and belt widths, with radius histograms. The first run (700 systems) is in the project reports: moons are 79 to 98 percent of the Class D bodies, nothing smaller than 215 km exists, no free-orbiting dwarf planet exists and a belt has no members. It measures only; no generated output changes.
