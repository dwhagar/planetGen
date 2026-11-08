**Headline:** The backup page no longer ticks a part that held nothing (the repair file of a site with no galaxy yet), so "Restore the ticked parts" works on those archives.

### Changed

- `GET /api/admin/backups` marks a part `empty` (its dump said `present` false) with the reason in `note`; the archive page shows it unticked and disabled with the reason.
- Restoring a whole archive leaves an empty part out and lists it in `skipped`, instead of failing with 409 on the part's "empty" refusal.
- The browser test of back up, compare, restore and import reads the Back up now message for every part the archive holds ("Backed up Settings, Sector repair file in ..."), which fixes its red on train 17.

generation-output: unchanged
