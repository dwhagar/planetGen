**Headline:** The galaxy skeleton step is no longer stopped as "stuck" while it waits for the generator battery's digest, and the wait now shows real counts: the battery's cases on the bar and the sampling tables built in the log.

### Changed

- OPS.184 (a bug from Boss's run on 8.0.4043, 5:36 PM EST on 2026-11-13): the plan's skeleton unit reports nine phases since OPS.170, and the phase "Waiting for the generator battery's digest" stood on one count while the digest process worked out the battery and the bright-star sampling tables of a code key it had not seen (minutes on a slow host), past the 300 s quiet limit twice. The run stopped with "unit 0 of stage Galaxy skeleton was stuck twice".
- The digest process now counts its cases and every sampling table it works out in two shared numbers (`generation/battery_progress.py`). The plan's wait reads them every 2 s: the bar moves a case at a time between phase 5 and 6, the log says `Generator battery digest, cases: 3 of 8 done (37%), 40 sampling tables built, ...` every progress interval, and the unit's heartbeat is told (`unitbeat.moving`) whenever the counts change, so a digest process that works never looks quiet and one that stops still does.

generation-output: unchanged
