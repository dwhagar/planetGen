**Headline:** A galaxy's species, civilization ages, technology levels, governments and colors no longer depend on how many workers filled it.

### Changed
- The population pass keys every draw on the homeworld's unique ID instead of a row id: a species' traits and age (`life:<id>`), its name (`species-name:<id>`), its technology indices (`tech:<id>`), and its polity's government and color (`government:<id>`, `polity-color:<id>`). Row ids come from per-worker id blocks, so the same galaxy used to draw differently at 1 and at 4 workers and in another fill order. Planets are scanned in unique-ID order. A polity color is now a drawn hue, so neighbouring polities no longer step round the color wheel.
