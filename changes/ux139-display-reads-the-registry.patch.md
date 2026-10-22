**Headline:** The unit ladders take every factor from the unit registry, and a `Quantity` can be shown on the site's ladders directly.

### Changed
- **Ladder rungs and companions are built from the registry** (UX.139): `ladders.rung("km", start=0.0, ...)` and
  `ladders.companion("atm", ...)` read the factor, offset and name from `unit_registry`, so a factor is stated once;
  the ladder keeps its display data (rung order, thresholds, asides, singulars, unit systems, and a symbol where it shows
  another, "lunar masses" for M☾). The rendered text and the generated `unitladders.js` are unchanged.
- **`format_measure(quantity)`** shows a `Quantity` on the ladder of its kind (`ladders.for_kind`) from its base value.
  The `format_*` names for a unit other than the ladder's base (km, AU, ly, pc, km/s, years, atm, g, solar masses,
  solar luminosities, particles per cm³, gauss) now share one helper over the registry's factors instead of eleven
  typed factors; each gives the same text as before.
- `unit_registry` gains the millibar and takes the Fahrenheit offset from `constants.CELSIUS_ZERO_K`; `to_base` and
  `from_base` leave a value alone when the unit has no offset (a `-0.0` stays `-0.0`).
