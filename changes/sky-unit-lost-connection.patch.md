**Headline:** A generator run no longer dies when MySQL drops a sky-tier connection, and shows the real error when it does.

### Fixed
- **A run stopped at stage 13 with `InterfaceError: (0, '')`.** A sky-tier unit lost its MySQL connection after about half an hour, and the rollback on the dead connection raised over the real error. A rollback on a connection that is already gone is now skipped, and a sky unit goes on from its last committed batch on a new connection (up to three).
