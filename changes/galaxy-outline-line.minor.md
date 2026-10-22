**Headline:** Planning a galaxy now says how dense it is and how much of it is planned, on the Generate page under the shape dials and in the plan output, and refuses a galaxy with nothing to plan.

### Changed
- **A line under the shape dials.** The Generate page shows, as the dials are typed, how dense the thin disc's crest is as a multiple of the one-star-per-sector density, about what share of the stars sit in planned sectors, and how many shells and candidate sectors the outline holds. A crest under 20 times the one-star density adds a warning.
- **The same line in the plan.** `planetgen plan` prints the same line, worked out by the same code, before it stores the outline.
- **An empty outline is refused with the reason.** A galaxy no sector of which reaches the one-star density stops the plan with what to change (raise the mass or the density, lower the scale length or the height), whether or not a dial was pinned.
