**Headline:** The stage-wording tests no longer expect "0 s" from a stage that really took a few milliseconds.

### Changed

- TEST.759: a stage line says how long the stage took on the units code's one duration ladder (UX.318), so a stage of 2.86 ms says "2.86 ms". The code is right; five tests in `test_stages.py` read the real clock and expected "0 s" or a whole-second "31 s". They now hold the clock the stages read still (`stage_clock`), and one new test pins that a stage of a few milliseconds says so.
