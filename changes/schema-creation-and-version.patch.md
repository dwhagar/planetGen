### Fixed
- **Several programs opening a new, empty database at the same moment no
  longer trip over each other (DB.5).** Each one used to create the tables
  itself and all but one could fail with "Duplicate entry ... for key
  'PRIMARY'". Now one creates them while the others wait, then carry on.
  The same goes for the control database (admin logins and the work
  queue), and for a migration running while another program connects.
- **A database whose schema version record was emptied or lost is no
  longer treated as up to date (DB.4).** Its version is now worked out
  from which tables and columns it has, so `migrateDb.py` (and update.sh)
  still runs the steps it is missing. No schema change.
