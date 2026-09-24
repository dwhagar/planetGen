### Fixed
- **Airless planets and moons could keep an atmosphere they no longer
  had.** A body moved into an airless class after changing zones kept its
  old class's atmosphere density, molar density and scale height. Those
  values are now cleared when it's reclassified.
- **Bodies around very dim stars could be colder than space itself.**
  Surface temperatures are now floored at the cosmic microwave background
  (2.725 K).
- Schema v30 cleans up both problems in rows already saved: stale
  atmosphere values on airless bodies go back to empty, and temperatures
  below 2.725 K are raised to it. Run `migrateDb.py` (update.sh does).
