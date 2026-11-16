**Headline:** A sector unit that runs again after its first try already saved the sector no longer fills it twice, which stopped a new-galaxy run on a duplicate object ID.

### Fixed

- DB.131: Boss's new galaxy (8.0.4317) stopped in the neighborhood stage on `IntegrityError 1062 Duplicate entry ... for key 'star_systems.uq_star_systems_uid'` (sector 75B000100861, serial 0, body 0). A sector is saved, and its cell charted, in one transaction, and the job engine marks the piece done only afterward; a unit run again after the save committed (a connection lost at the commit, a worker killed just after it) filled the cell a second time, and the second fill's systems took the object IDs the first fill's already held. `run_galaxy.generate_sectors_do` now asks the chart ledger first: a charted cell is reported from the rows (`store.saved_sector_report`: sector id, name, systems, stars, planets, phenomena, `stored` true) and is not filled. The sectors it reports add to the run's counts and are left out of the speed figures (`run_common._count_reported`). A held run that stopped on this error is finished by Retry held once the fix is installed. Reuse: `chart_ledger.row` says whether a cell is charted, as `generate_sectors_pieces` already does for a run begun again.

generation-output: unchanged
generation-affects: none
