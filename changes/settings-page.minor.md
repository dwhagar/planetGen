**Headline:** Admins can now see and change every program option from a new Settings page, with help for each option, a review of the changes before anything is saved, a backup of the old file and an audit row for every change.

### Changed
- New Admin page `/admin/settings` (a Settings tab on every admin page and a link on the Admin overview). It lists every option of `config.json` grouped by category, each with its help text, default, limits, unit, where the value in use comes from (default, `config.json` or an environment variable) and whether a change needs a restart. The list is drawn from the settings model, so an option added to the model appears with no page code.
- Press Review the changes to see the old and new value of each changed option, which ones need a restart and which an environment variable still overrides; a value that does not fit is sent back with its message at its option and nothing is saved. Confirming saves through the settings writer: a backup archive of the old file first, an atomic write and one audit row per option. A file changed by someone else after the form was opened is refused, not overwritten.
- Options an admin may not change (command line only, or the Owner's) are listed read-only with the reason. A secret is never shown: its box is write-only. Use the default removes a key from the file.
- While a saved option still needs a restart of the web server, every admin page shows a banner naming it.
- New API routes `GET /api/admin/settings`, `POST /api/admin/settings/plan` and `POST /api/admin/settings` (`docs/api.md`).
