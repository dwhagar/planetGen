### Changed

- Stars and phenomena draw their generated name only when it is first read (PERF.43). A placed phenomenon is named by its object ID and never draws one, which saves about 4 percent of a sector fill; building a rogue planet is about four times faster. The same seed gives a different galaxy once, since each name now takes one draw at construction instead of many.
