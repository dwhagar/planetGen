**Headline:** A browser test's database can no longer hang the test run for fifteen minutes when it is torn down.

### Changed
- **Test servers stop cleanly (TEST.193).** `tests/served.py` runs the web server of a browser test file and, when the file ends, waits up to 8 seconds for requests a page had already sent, so none still holds a database connection when the test database is dropped. The Galaxy Map's tile requests held one for 4 to 6 seconds under load.
- **Dropping a test database is bounded (TEST.193).** `dbpool.drop_database` waits at most 10 seconds for a lock, ends the connections that hold the database, and names each in a warning with the code that took it, so a leak shows up in the run instead of as a 853 s wait. Test-only; nothing in the program changed.
- **One fewer fetch under load (TEST.195).** The sector browser tests work out the sector's map designation in the test process instead of fetching the "Show on Galaxy Map" redirect, a request that timed out at 30 s while the page's own tile requests were running. The redirect itself is still tested in `test_web_sector_nav.py`.

TODO items resolved: TEST.193 (cause found, fixed and bounded), TEST.195 (the fetch removed; the 30 s stall itself was never seen again, see the issue). Added: none.
