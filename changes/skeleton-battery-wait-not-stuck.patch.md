**Headline:** The galaxy skeleton step is no longer stopped as "stuck" while it waits for the generator battery's digest, which takes minutes the first time after an update.

### Changed

- The plan's skeleton unit reports its nine phases since OPS.170, so its counter can be judged; the phase "Waiting for the generator battery's digest" stands on one count for as long as the digest takes in its own process (minutes when the digest and the bright-star tables of a new code key are not kept yet), which passed the 300 s quiet limit twice and stopped the run ("unit 0 of stage Galaxy skeleton was stuck twice"). The wait, and the end of that process, are now announced to the heartbeat as a wait on the generator battery, which has its own 30 minute flag and is never judged quiet.

generation-output: unchanged
