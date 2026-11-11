**Headline:** The Queue pages show run durations, ETA ranges and heartbeat ages in their original form again ("3 m 05 s", "12.4 s").

### Fixed
- **OPS.142 defined a second `_duration` in `web/queue_page.py`**, which replaced the first at import and changed every figure that used it. The other-job rows' short form is now `_beat_gap`; a test checks the two forms and that no function of the module is defined twice.
generation-output: none
