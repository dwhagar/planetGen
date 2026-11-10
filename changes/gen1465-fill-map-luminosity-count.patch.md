**Headline:** `planetgen fill-map-luminosity` no longer reads the whole scatter table once before its first unit starts.

### Changed

- The expansion of the fill counted `phenomenon_scatter` rows with a mass and no map luminosity yet to print how many there were. No index serves that condition, so on a full galaxy the count was a read of the whole table, single-threaded and outside the memory gate, before any unit ran. It is gone; the plan line names the highest id instead. The units read the same rows, a bounded id range each, in parallel and under the gate.

generation-output: unchanged
generation-affects: none
