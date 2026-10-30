**Headline:** The map-luminosity fill and the habitability re-score now run as stages of the shared run, with their finishing steps inside it.

### Changed
- **The two fills no longer start a run of their own (PERF.202).** `planetgen fill-map-luminosity` is a map-luminosity stage followed by the Galaxy Map's block light totals made again, as the scatters already were; `planetgen rescore-habitability` is a habitability stage whose block totals are made by the command's own finishing step, once the last range is done. The output is unchanged, and the stages can now share one run with other commands.
