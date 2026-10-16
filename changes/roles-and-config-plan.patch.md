**Headline:** Version 10.0 is planned as roles and permissions (user, moderator, admin edit requests), and the configuration page moves up into 8.1.

### Changed

- The configuration page (ADM.43, ADM.44) moves from 8.2 to 8.1 phase 1 as its own thread, broken into ADM.83 to ADM.85; its Owner tier is ADM.86 in 8.2.
- 8.6 plan items (DB.36, ADM.58, ADM.60, ADM.66, ADM.68) carry the requester, the approving role and the `requested` state, so 10.0 needs no new plan table; an admin confirms their own plan with no second approver.
- USR.2 stores four roles (user, moderator, admin, Owner) behind one guard so the roles of 10.0 fit without a rebuild.

### Added

- The 10.0 headline USR.25 and its items USR.12 to USR.24, DB.39, TEST.178, TEST.179 and DOC.55, the design note `docs/design/roles-and-requests.md` and the version file `docs/plan/version-10-0-roles-and-permissions.md` in four phases.
- The per-cell system cap as a configuration option (GEN.303), wired to the settings page.

TODO items resolved: none / Added: ADM.83, ADM.84, ADM.85, ADM.86, GEN.303, USR.12 to USR.25, DB.39, TEST.178, TEST.179, DOC.55
