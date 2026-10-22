**Headline:** A Generate run no longer sits for hours on one unit drawing a Magellanic Cloud: a large satellite or cluster is drawn in slices that the workers share, and every unit of components reports its progress.

### Fixed
- **A Generate run stalled for hours on its components (the star scatter's last unit, then the phenomena's).** All the globular clusters, open clusters and satellite galaxies were one unit, drawn by one worker one star after another. With a low luminosity floor the Large Magellanic Cloud alone holds about a million stars, so that unit was nearly all of the stage's work, took more than an hour on a single worker while the other workers stood idle, and reported nothing until it ended, so the bar stuck near the start with an estimate of 30 hours. The components are now cut into tasks by their expected stars (or phenomena) alone, the same at any worker count: a component with more than 25,000 expected stars (1,000,000 expected phenomena) is cut into slices, one unit each, and the smaller ones are bundled in order up to that size. The units of components report the stars they have drawn a few times a second, so the bar and its estimate follow them.

### Changed
- **A sliced component is drawn from streams of its own.** A slice draws from the component's stream with `:slice:` and its number added, with its share (1 over the slices) of the component's Poisson mean, so the slices together are a draw of the whole component and a seed gives the same stars at any worker count. A component that needs no slices keeps its stream and its stars; the Large and Small Magellanic Clouds, the Sagittarius dwarf and any other component past the limit at the galaxy's floor draw different stars from before.
- **The plan of a scatter names its component tasks.** `component_key` of the shell scatter's and the phenomena's plans became `component_parts`, task by task the `[kind, name, slice, slices]` lists of the components it draws.
- **A scatter's log lines say which unit they are about and how big it is.** A line now reads `Bright stars, shells 100-131 (unit 12 of 512, about 1,200 stars expected): added ...`, and the components' and the phenomena's the same, so the sizes of the units can be read off the log.

generation-output: changed

TODO items resolved: PERF.148. Added: none.
