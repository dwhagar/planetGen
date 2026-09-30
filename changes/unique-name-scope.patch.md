### Changed

- Renames check only uniquely named objects (sectors, systems and stars)
  for a clashing name. Planet and moon names come from their star's, so
  `PATCH /api/planets|moons|stars|systems/<id>` no longer searches the
  planet and moon tables; a planet or moon may share another body's name.
  Generation already skipped them (since v34). Profiling a galaxy run
  showed name reservation at about 1% of generation time; row-by-row
  inserts are the main cost.
