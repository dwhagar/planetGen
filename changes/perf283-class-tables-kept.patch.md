**Headline:** The site's class catalog no longer rebuilds its two star tables in every process, so the engine's app template starts about 13 s faster on a cold box (PERF.283).

### Changed

- The forced-star tables (one per [Fe/H] bucket) behind the class reference catalog are kept in the sampling-table cache the bright-star tables use, so the app template, the web workers and the warm-up build each once. Measured: the class catalog took 14.3 s the first time and 1.4 s read back, with the same cells.
