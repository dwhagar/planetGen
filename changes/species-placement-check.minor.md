**Headline:** A species is only placed where it could live, and a report lists any that are not.

### Added
- **The placement check (POP.14).** `species_fit.check_placement` is the one check every way of putting a species on a planet or moon calls before it writes. It refuses a placement outside the species' limits and says which conditions fail in plain words ("it is too hot"). A habitat or vessel may be placed where the species could live with the equipment its tech level can make. A species on its own homeworld is a native there, so its Suitability is not judged.
- `python3 -m planetgen.cli.species_check` lists every species placed outside its own limits (exit 1), or says all are inside (exit 0).

### Changed
- The population pass makes no species on a homeworld whose own limits cannot hold it (undrinkable water); that planet is passed over like any planet without a civilization.
generation-output: changed
