**Headline:** A supernova remnant on the Sector Map is now drawn from its age: young ones a thin bright rim, old ones a thick faint shell.

### Changed
- **Remnant shells by age (MAP.286).** The sector scene now carries each supernova remnant's age, and the shell is painted from it: a 50-year-old remnant is a thin, bright, hot-coloured rim around a hollow centre, and a 100,000-year-old one a thick, faint, cooler shell with the hollow filling in, with every age between interpolated on a log scale. The remnant's card shows its age, and the sector API's `phenomena` entries gain `age_years` (null for other types). The look is a stand-in that follows the body-rendering design table; the filament texture and colour by emission class are still to come. Nothing generated or stored changes.
