**Headline:** The reserved migration revision numbers are recorded as an order of merge, not reservations.

### Changed

Core 3 found that the migration runner requires contiguous revision ids, so a reserved revision number cannot be used until every lower number exists. The plan notes on GEN.324 and DB.53 now say that the listed numbers are an order of merge and each migration takes the next free number when it merges.

TODO items resolved: none / Added: none
