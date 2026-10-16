**Headline:** OPS.51 is retired (the maintenance switch now lives where the web user can write it), and a Liquid delimiter in a TODO note is removed.

### Changed

- OPS.51 (#1741) leaves the open list: the maintenance switch file lives in /var/lib/planetGen/state and a switch that cannot be written fails loudly.
- TEST.163's note no longer holds a Liquid delimiter, so `test_docs_liquid` passes on main again.
- Interface is idle and its queue is unassigned in the execution plan.

TODO items resolved: OPS.51 / Added: none
