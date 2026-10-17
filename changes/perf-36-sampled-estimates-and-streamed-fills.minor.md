**Headline:** A fill of tens of thousands of sectors or more is estimated from a sample and generated from a stream, and a web request is refused a ball of more than 2 million sectors.

### Changed
- **Large estimates come from a sample (PERF.36).** A run of more than 20,000 sectors (a shell, a span, a course or a radial fill) works its size and time out from 200 sectors spread evenly through it and scales them up; the line says "worked out from 200 sampled sectors". Smaller runs are worked out from every sector as before.
- **A large fill is no longer held in memory (PERF.36).** Its sectors are counted in a first pass and resolved one by one as the work queue takes them, instead of being listed first.
- **A web request cannot list a huge ball (PERF.36).** The sector page's neighborhood action answers `413`, and the Generate page's radial estimate stops, when the radius holds more than 2 million sectors; the command line is not limited.
