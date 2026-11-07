**Headline:** The update guard no longer believes what a local user leaves in the shared run and stop-reason folders.

### Changed

- `admin/updateguard.py` (SEC.38, audit finding S-M2): a marker in `running/` is listed and canceled only when it is a plain file named `<kind>-<pid>.json` for its own pid and owned by the user the process of that pid runs as; a marker that fails this, a link and a pipe are ignored and never opened, so a local user can no longer have the guard (root, in `update.sh`) signal another user's process. A run registers its marker on a name that does not exist yet (`O_EXCL`, `O_NOFOLLOW`, the name cleared first).
- `admin/stop_reason.py`: a reason is written the same way, so a link planted under `<pid>.txt` is replaced and never written through; a run reads a reason only from a plain file owned by root or by itself, never through a link.
- Tests for the forged marker, the mismatched name, the link, the pipe, and the planted link on both write paths. The folders stay mode 1777 (every user a run starts as must be able to register); nothing in them is trusted.

generation-output: unchanged
