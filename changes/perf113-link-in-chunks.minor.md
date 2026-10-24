**Headline:** The neighbour link of a large run is cut into units of 500 sectors, so no single unit runs for long.

### Changed
- **The link stage of a run is cut into units of 500 sectors (PERF.113).** A run linked every sector it made to its neighbours in one unit (70 ms a sector at the default galaxy size, so one unit of several minutes for a run of thousands of sectors). The sectors are found once, put in neighbourhood order and linked in chunks. The links stored are the same; the stage is not faster, because linking takes the neighbour lock one batch at a time.
