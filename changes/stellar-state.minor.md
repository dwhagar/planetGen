**Headline:** The generator can now name a star's evolutionary state and the age since a white dwarf stopped burning.

### Added
- **Evolutionary state and age since it stopped burning (GEN.785 and GEN.786, definitions).** `physics/stellar_state.py` maps a star's Yerkes class and phase to one of twelve state words (main-sequence, red-giant, helium-burning, white-dwarf and so on) and gives a white dwarf's age since it stopped burning. Nothing in generation calls it yet; the wiring rides the star epoch.
