**Headline:** An update whose pull changed update.sh now runs the new script, not the old one it had already opened.

### Fixed
- OPS.301: bash keeps running the copy of `update.sh` it opened, so a pull that changed the script left the whole update on the old one (Boss's update at 6:05 PM EST 11-23 ran the math check that OPS.286 had taken out, and printed "of 9" stages for an 8-stage script). Right after the pull, `update.sh` now (checksum by `git hash-object`, which macOS has too) runs itself again from its start when it changed, under the same process id so the update lock stays his, and tells the second run the guard the first took (nothing is asked or reset twice). The update that first brings this still runs the old script; the one after it is the first to use it.

generation-output: unchanged
