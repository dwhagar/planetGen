**Headline:** A warm-up that a change started no longer runs while a generation run is live.

### Changed
- The warm-up a data change starts by itself (`planetgen warm --on-change`) now waits for any live generation run to end before it fetches pages, so it no longer logs failed warm-up requests against data the run is rewriting. It keeps its lock fresh while it waits, and the generation run's own end asks for the warm-up that follows. Output of generation is unchanged.
