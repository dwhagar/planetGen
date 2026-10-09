### Changed
- A `galaxy` run now links its new sectors to their neighbours (containment in nebulae and remnants, nearest systems) once at its end instead of one sector at a time while it saves, so the workers no longer queue for it (PERF.45). The stored links are the same. Sectors made on demand still link as they are saved.
- A planet's or moon's position is worked out once when first read instead of after every move (PERF.46); the values are unchanged.
