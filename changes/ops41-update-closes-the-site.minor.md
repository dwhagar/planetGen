**Headline:** The site is closed for visitors (maintenance mode) while an update runs, and opened again when it has finished.

### Changed
- **`update.sh` closes the site before it changes anything.** It turns maintenance mode on (cause "update") and takes a lock, so visitors see the "being updated" page instead of pages that hang while a migration holds the database. A good update opens the site again at its end. If you had closed the site yourself, it stays closed afterwards, and a sector-grid conversion keeps it closed until the new Galaxy Map is built, as before.
- **A failed update leaves the site closed and says why.** The message names the step that failed and the command that opens the site again (`python -m planetgen.cli.maintenance off`), because the database may be half migrated. The lock is released so `./update.sh` can simply be run again.
- **A migration run by hand is guarded the same way.** `python -m planetgen.cli.migrate` closes the site while it has something to migrate and opens it when done; with nothing to migrate it changes nothing.

TODO items resolved: OPS.41 / Added: none
