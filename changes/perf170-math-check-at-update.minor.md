**Headline:** Generation no longer spends 2.5 seconds on the math check before every run; the check runs at every update and from the Admin page.

### Changed
- **The math check no longer runs before generation (PERF.170, Boss).** `galaxy`, `plan`, `population`, `sector --num-sectors`, the Generate page's jobs (their "Check the math" step is gone) and the Sector page's "generate the neighbourhood" button start without it. `update.sh` still runs `planetgen check-math` at every update (and warns, and skips the population pass, when it fails), the website still runs it once at start-up, and `planetgen check-math` still works by hand.
- **The Admin page has a "Check the math" button.** It runs the checks when an admin asks and shows the summary, naming the failed checks.
