**Headline:** The admin Update page says plainly when the last update failed, and the update log and a queue run's log are colored by level (errors red, warnings amber, info blue, debug green).

### Changed

- `src/planetgen/admin/webupdate.py` (UX.284): when no update is waiting or running, the view now has `outcome` (`failed`, `done` or none) and `outcome_detail`. It is `failed` when systemd holds `planetgen-update.service` in `failed` (read with an unprivileged `systemctl show`: exit status 1 from the update guard, or the start limit of ten tries in ten minutes, which needs `systemctl reset-failed planetgen-update.service`), or when the last run in `update.log` ends in an error line. It is `done` when the last run reached its closing "Update done" line. A start that dies before `update.log` begins now shows too.
- `src/planetgen/web/templates/admin_update.html`: a red "The last update failed." alert under State, with the detail and the last error line of the log.
- `src/html/static/style.css`, `partials/log_lines.html`, `logrecord.display_severity`: each line of the Update page log and the queue run log is in its level's color (error and critical red, warning and notice amber, info blue, debug green), 4.5:1 or better on the box in the light and the dark theme. The level word stays in every line, so color is never the only cue.
- `docs/reference/api.md`: `GET /api/admin/update` also returns `outcome` and `outcome_detail`.

generation-output: unchanged
