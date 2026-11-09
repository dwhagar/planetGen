**Headline:** A running unit's state (working, slow, waiting, quiet, looping or dead) is now one tested rule set, and a unit whose counter is moving is never judged for killing.

### Added
- **The heartbeat judge (OPS.137).** `queue/beatjudge.py` turns a unit's beat record, its command's recorded gaps and a clock into a state with the rules of the heartbeat research: dead after 30 s without a beat, waiting on the gate or the sky lock, quiet when the counter stands still past `max(5 min, 5 x the longest recorded gap)`, looping when the count runs past 3 or 10 times its estimate, goes backward or restarts a phase, and slow below a third of the expected rate. Only dead and quiet may be acted on. Nothing calls it yet; the beat record (OPS.136) and the follower (OPS.138) come next.

generation-affects: none
