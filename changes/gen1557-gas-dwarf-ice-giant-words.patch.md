**Headline:** the gas dwarf and ice giant class texts now describe what the bodies are made of (rock and metal under a thick envelope), and the no-surface line says the rock lies deep under it.

### Changed

- Classes T and I in `physics/data/planet_classes.json`: the gas dwarf was "a deep hydrogen and helium envelope over a small core" and the ice giant "mostly water, ammonia, and methane ices"; the stored bodies are rock and metal 80 to 90 percent of the mass under an envelope of 1 to 18 percent (class T), and two sampled ice giants hold no water. The gas dwarf now reads "mostly rock and metal under a thick hydrogen and helium envelope, with no surface to stand on"; the ice giant "a thick hydrogen and helium envelope over an interior of rock, metal, and the ices its disk gave it". The class page, the planet page and the system map read these words when they are built, so nothing stored changes.
- The day-counter line of a body with no surface (state word `no solid surface`, unchanged) now reads "No surface to stand on: any rock lies deep under a thick envelope. The day below counts the turn of its clouds."
- New test `test_giant_class_words.py`: the texts name rock and do not say "small core" or "mostly ices", and the gas dwarfs of a generated system store more rock and metal than envelope.

generation-output: unchanged
generation-affects: none
