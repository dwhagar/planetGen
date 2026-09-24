### Changed
- **System pages are rendered natively from the database.** The System
  panel is now an expandable list of the system's stars, planets (moons
  nested under each), asteroid belts and comets, each row showing class,
  terrestrial/gas giant, Habitable yes/no and Inhabited yes/no, and
  opening onto that body's own description. Wikitext and Markdown buttons
  show the full generated page in a code box with a Copy button.
- **Wiki page text is no longer stored (schema v29).**
  `star_systems.wikitext_content`/`markdown_content` are dropped; both
  formats are rendered on demand from the system's rows
  (`stellarObjects/systemRender.py`), so pages now follow renames, names
  made unique after generation, and orbit ticks. Wiki upload uses the
  fresh render. `GET /api/systems/<id>` no longer returns the two text
  fields; use the new `GET /api/systems/<id>/text?format=wikitext|markdown`
  and `GET /api/systems/<id>/sections`. **The migration deletes the stored
  copies: back up first** (`mysqldump`, or `src/checkRenderParity.py
  --export-dir`, which also compares the stored and rendered text on a
  database still at v28).

### Fixed
- **Binary stars could load with their primary and secondary swapped**
  when the second star generated was the heavier one, changing the
  "Barycenter Offset", planetary-limit and "Location" lines on reload.
