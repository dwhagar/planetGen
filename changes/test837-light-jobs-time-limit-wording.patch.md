**Headline:** The light job tests expect the time-limit message in the wording the product writes (1m 40s, 1s), so they pass again.

### Fixed

- TEST.837: `test_light_jobs.py` expected the time limit of a job that ran too long as 100 s and 1 s; the product writes it as 1m 40s and 1s since the duration rule. Tests only.
