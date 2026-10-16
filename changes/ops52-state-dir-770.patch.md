**Headline:** The state folder stays group-writable (770) however the permissions script is run, matching create-cache-dir.sh.

### Fixed
- `set-permissions.sh` leaves `/var/lib/planetGen/state` at 770 (owner and group may write, others nothing) instead of 750, so re-running it alone no longer takes the group's write bit away from an operator running `python3 -m planetgen.cli.maintenance`.

TODO items resolved: none / Added: none
