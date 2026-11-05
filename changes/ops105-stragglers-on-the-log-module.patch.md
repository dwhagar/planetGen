**Headline:** The activity log, the queue run's live log and the last direct users of Python's logging now go through the one logging module.

### Changed

- `util/log.py` makes the activity log's file handler (with its file mode and its own rotation), the run log's handler and the sink the activity lines go to; `admin/activity_log.py` and `queue/runrecord.py` no longer define a handler, and four modules that fetched a logger by name now call `log`. The activity log's lines look as they did after the local-time change (`2026-10-01T04:00:00-0400 planetgen[1234]: AUTH login.failed ...`), so the fail2ban filter is unchanged by this.
- The test that watches for stray formatting now also fails on any handler, formatter or `logging.getLogger` outside `util/log.py`; its list of allowed strays is empty.
