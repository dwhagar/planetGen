**Headline:** Two test fixes: the runaway-velocity test no longer fails now and then, and the plan tests draw a naming key instead of logging that they cannot.

### Changed
- **The runaway-velocity test picks a system with a body that moves (TEST.191 and a flaky test).** Its setup stopped at the first system with any planet, but a system whose only entry is an asteroid belt has no body with a velocity, so the test failed about once in 300 runs. It now waits for a system with a body that is not a belt.
- **The plan tests have a control database to draw the naming key into (TEST.191).** The real-plan tests no longer log "No naming key drawn: the control database can't be opened"; the control schema is set up in the test database, and a test checks the key is stored. The other half of the item, stray files in the jobs directory's `galaxy-settings` folder, went with ADM.50, which removed that folder.
