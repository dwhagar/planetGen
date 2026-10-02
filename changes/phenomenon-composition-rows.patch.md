### Fixed
- **Asteroid field and interstellar comet pages now show their saved
  composition (DB.2).** Each component was saved in its own row but never
  read back; the page showed the one-line summary stored beside it. The
  page now builds the Composition line from those rows, and
  `GET /api/phenomena/asteroid_field/<id>` and
  `/api/phenomena/interstellar_comet/<id>` return them as `composition`
  (a field's as `{component, concentration}`, a comet's as a list of
  components). No schema change.
