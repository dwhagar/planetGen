**Headline:** The update's warm-up progress bar follows the work inside each unit, and a keypress skips the warm-up.

### Changed
- **OPS.106: update.sh no longer guesses how long a step takes.** The "a few seconds", "5 seconds to a minute" and "a few minutes on a big galaxy" texts are gone; the progress bar says how long is left.
- **OPS.107: the warm-up bar moves inside a unit.** Sector, system and phenomenon slices report how many objects they have warmed, and the bar, the stage's time left and the whole warm-up's time left count that part of a chunk, so a big galaxy's few large chunks no longer leave the bar at zero until the first one ends. `planetgen warm --follow` shows the same.

### Added
- **OPS.108: press s during the warm-up to skip it.** update.sh shows the key while a warm-up stage runs (only with a terminal; without one nothing is read and the warm-up runs as before). What was stored is kept and the first visit builds the rest.
- **OPS.109: the data set download step lists what it keeps by download URL.** Files already stored are shown as one `kept` line for each URL (with the file count), not one line for each file, and the dataset line says the folder they are stored in.
