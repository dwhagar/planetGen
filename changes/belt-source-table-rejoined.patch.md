**Headline:** the belt-source stand-in table is whole again, so its fourteen lost rows are counted and checked

### Changed

A merge had left the belt-source table (6.95 of the gap-filling plan) with a duplicate of an old BR8 row, an old BR9 row that the roundness work retired, and a blank line between them and BR10 to BR23. The registry reader stopped at the blank line, so those fourteen rows were missing from the stand-in values list and one test (the burn-down register, row BR20 at L3) saw its line as stale. The blank line and the two old rows are gone, the list is rebuilt (1,044 rows), and the two body-status rows that now show (BR22 and BR23) have wiring lines under GEN.1092. No generated output changes.
