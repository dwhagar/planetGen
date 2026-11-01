**Headline:** every stage with units now holds its workers to what a worker and its horse really take of the host's memory

### Changed

A worker and its horse measured 2.4 GB for a sector and 1.5 to 1.6 GB for a link, a settle or a map-luminosity fill (the horse alone showed 1.5 GB for a sector, because its resident size counts the worker's tables but not the worker beside it). The largest-worker figure in the run report is now the two together, shared pages counted once, and the figures are 2,600 MB for the sector stages and 1,700 MB for every other stage with units, the scatters included, so a 16 GB host runs 3 sector workers or 4 of the others at once where it ran 36 and was killed. A worker freezes what it holds before it forks a horse, so the garbage collector does not write to every page of its warm tables in each horse. No generated output changes.
