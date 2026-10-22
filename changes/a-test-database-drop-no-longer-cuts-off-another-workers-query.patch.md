**Headline:** Dropping a test database no longer ends another worker's query that is counting its rows, which was the "Lost connection to MySQL server during query" failure under load.

### Fixed
- **Tests failed with "Lost connection to MySQL server during query" or a 503 "database unavailable" while four workers ran (TEST.262, the sightings listed on it).** The drop helper of the test suite ended every connection holding a lock on a table of the database it dropped, and another worker's `list_databases` (the database picker and the admin pages) holds exactly such a lock for a moment while it counts the rows of every test database. The helper now ends only the connections whose current database it is, gives any other connection three seconds to let go, and ends it only if it does not. A loop of three listers beside three create-and-drop loops failed 222 times in a minute before and 0 times after.

TODO items resolved: none. Added: none.
