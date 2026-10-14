### Fixed
- A fresh install no longer breaks every page that shows rendered Markdown (system descriptions, the population text): `markdown` 3.11.1 removed two inline patterns the renderer named as required, and the renderer now removes only the ones the installed release has.
