**Headline:** The memory controller now learns where its pressure edges belong from what it read in earlier runs: it records each minute, and moves its edges up or down when memory kills came unseen or holds were never needed.

### Changed
- **The controller's edges are learned from a record (OPS.249, H4 of the job engine; answers OPS.151's question).** The supervisor writes one row a minute of what the memory controller read (`pressure_minutes`, control schema v25: the highest `some` and `full`, the lowest available memory, the highest level, the seconds admission was held, the jobs asked to stop). `queue/pressurelearn.py` reads the last 14 days and the ledger's `killed_memory` attempts and works out one sensitivity, from half to double the typed one: the PSI edges are divided by it and the shares of the reserve multiplied by it, so the ratios between the edges stay as they are. Three or more memory kills with no pressed reading before them lower the edges; eight or more holds of which 80 percent ended with no starved reading, no job asked to stop and no kill raise them, by up to half. With no record, too little of it, or `PLANETGEN_PRESSURE_LEARN=0`, the typed edges stand. A move is said in the supervisor's list.
- The database signal's factors and floors (`dbsignal`, OPS.248) are not learned in this change.

generation-affects: none
