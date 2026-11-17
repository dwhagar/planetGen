**Headline:** The supervisor now holds the bulk and upkeep lanes when the database's own answers slow down, as it already does when memory runs short.

### Added
- `queue/dbsignal.py`: a thread times one ledger read every two seconds; the calm time is learned from the database itself; a read at four times calm (and over 50 ms) presses, ten times (and over 250 ms) starves; a read still out counts as the time it has been out.
- The hold is the memory controller's own machine with this signal for its look (reused, not copied): after ten seconds in a row nothing of the bulk or upkeep lanes is admitted and thirty seconds of calm lets it go, each change said in the supervisor's list.

### Changed
- `pressure.Controller` takes the reason it gives for a hold; the memory controller's words are unchanged.

generation-output: unchanged
