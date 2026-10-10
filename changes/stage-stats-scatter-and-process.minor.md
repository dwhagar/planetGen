### Changed
- Every staged process now stores its whole-process time, even a run that visited no layers, and each stage's seconds per layer. The scatters are stored on their own (massive stars, bright stars, phenomena, backfill) and added up in one "All scatters" row.
- A scatter with no stored time of its own is now estimated from that common scatter figure, in the command-line bar, the Generate page's whole-job bar and each stage's bar alike, since they all read the one stored-time function.
