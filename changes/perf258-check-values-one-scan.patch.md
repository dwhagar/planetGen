**Headline:** The database check for impossible values reads each table once instead of once for every rule, 34 scans instead of 5, with the same findings.

### Changed

- **PERF.258: `check-values` counts all the rules of a table in one scan.** It had 34 rules (13 on planets, 13 on moons, 5 on stars, 2 on star systems, 1 on sectors) and ran a separate count over the whole table for each, so planets and moons were each scanned 13 times and stars 5 times. It now makes one pass per table with a conditional sum for each rule, and reads a rule's rows again only when the rule finds some, to list them. The problems, their order, their wording and the sector scope (`--shell`, `--sector`) are the same, and a test compares every rule's count with a count of that rule alone and checks that exactly five scans are made. No generated output changes.

generation-output: unchanged
generation-affects: none
