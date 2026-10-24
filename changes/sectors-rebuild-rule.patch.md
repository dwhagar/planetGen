**Headline:** The rule for changing the `sectors` table is written down: batch the changes, say the rebuild cost, keep the name index.

### Changed
- **Decided how later schema changes pay for the `sectors` FULLTEXT index (DB.84).** The index forces a full table rebuild for any column change, so the migration helpers' docstring, `migrations/README.md` and the architecture note now say: keep the index, batch every `sectors` change of a version into one revision, and say in that revision's docstring that it rebuilds the table. No code behaviour changes.
