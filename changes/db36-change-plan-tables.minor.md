**Headline:** The database gets the tables the 8.6 edit-everything plans will use: plans, their operations and the edit journal.

### Changed
- **`change_plans`, `change_plan_ops` and `edit_journal` (DB.36, schema v122).** A plan records who asked and in which role, the request, its state (draft to applied, undone, stale, rejected or withdrawn), the role that must approve it, who decided and when, and the job that applied it; its operations are rows in order with the target object ID, the values before and after as JSON and a status; the journal is the append-only list of every applied operation with the object, sector, aspect, before, after and generator epoch, keyed by object and by sector. The requester, approver and state fields are there from the start so the 10.0 roles need no rebuild. The migration only makes the tables; nothing reads them yet.

TODO items resolved / Added: DB.36 / Added: none
