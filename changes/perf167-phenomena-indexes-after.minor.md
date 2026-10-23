**Headline:** The phenomenon scatter of a large galaxy finishes in 73 minutes instead of running past 90, because two more indexes wait until its rows are in.

### Changed
- The phenomenon scatter now leaves the position index and the due-time index off while its rows are written and builds them after the last row, as it already did for the address and class indexes. At the large size (63 million phenomena) the step took over 90 minutes (stopped at the cap) before and 72.7 minutes after, with half the InnoDB writes. The rows are the same.
- `docs/design/generation-performance.md` has the measurements and what is left.
