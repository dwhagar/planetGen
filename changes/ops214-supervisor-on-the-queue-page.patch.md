**Headline:** The Queue page now shows whether the supervisor is alive, since when, and what it did lately with the reason.

### Changed
- **OPS.214, output unchanged.** A panel "The supervisor" on the admin Queue page (`GET /api/admin/work/supervisor`): alive or not, since when, how long it has been up, how long since its last turn, its turns and process, and its latest twenty acts newest first, each with its time and its reason in words (`started: something asked for one`, `ended: nothing to watch for 300 s, after 150 turns`). `Supervisor.act(what, reason)` is the one call every later item uses to say what it did, which also goes to the log.
generation-output: unchanged
generation-affects: none
