### Added
- **`docs/server-checklist.md`**, a step-by-step check to run on the server
  after a deploy: confirms the code version, that `migrateDb.py` has brought
  every database to schema v26 (an Apache restart alone doesn't), the
  spatial indexes, `request-timeout=60`, the Galaxy Map tile cache, and that
  the Galaxy Map no longer times out or runs the API out of memory.
