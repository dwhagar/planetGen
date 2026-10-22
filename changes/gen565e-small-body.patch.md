**Headline:** Comets, interstellar comets and asteroid fields now share one `SmallBody` family class.

### Changed
- **`generation/small_body.py`** (GEN.565e): `SmallBody(LazySaladName, MassObject)` holds the constructor preamble the
  three small bodies each wrote out (the shared configuration, a name drawn on first read, whether the name was given) and
  the save and load of a composition kept as a list of strings. `Comet`, `InterstellarComet` and `AsteroidField`
  subclass it; the field keeps its own pair-list save and load. No field, stored value or draw changes.
