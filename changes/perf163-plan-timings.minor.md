**Headline:** The plan prints where its time went, and the plan of a default-size galaxy takes about a fifth of what it did.

### Added
- **The plan log ends with a table of where the plan's time went**, costliest step first: the galaxy's shape, the globular clusters, the satellites, the young open clusters, the component field, the outline, storing the components, the shells and the shape.

### Changed
- **The plan no longer works out the clusters' system counts twice.** The cache of the mean system mass at an age held 512 ages while a plan asks for about 22,000, so the line that prints the counts paid the whole cost again; it now holds 65,536. The search for the mass of a given lifetime also reads the metallicity's terms once. No drawn value changes.
