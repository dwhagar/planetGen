**Headline:** The one-off system page now takes the number of stars, says that habitable world and intelligent life only raise the odds, and refuses bad input with the same messages as the command line.

### Changed
- `/admin/generate/system` has "Number of stars" and "Most stars allowed" boxes (`--stars`, `--max-stars`), checked against the pair options and the system file the way `planetgen system` checks them.
- Habitable world and intelligent life are offered as "Raise the odds" instead of "Yes", since the generator never forces them (GEN.1130).
- The form now refuses, before running the generator, a star type that is not a spectral type, a name longer than the database allows, and no planets combined with an asteroid belt.

generation-output: unchanged
