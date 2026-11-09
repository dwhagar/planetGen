# Items: SEC, Security

The life of every SEC item that a commit title names, from the git history as of `d02661568`: first and last date, whether it is open in [docs/TODO.md](../TODO.md) now, and its commit titles in order (the pull request numbers are the GitHub pull requests). An item not open in TODO.md has shipped, been closed as a finding or been dropped; [dropped-and-closed.md](dropped-and-closed.md) lists the closings.

## SEC.1

2026-10-01 to 2026-10-01; 1 commits; not open.

- 2026-10-01: Lock out an address after failed logins (SEC.1, SEC.21) (#220)

## SEC.20

2026-10-01 to 2026-10-01; 2 commits; not open.

- 2026-10-01: TODO: bugs-first, then security plan; login brute-force research (SEC.20-SEC.27)
- 2026-10-01: Always-on activity log and failed-login records (SEC.28, SEC.20) (#217)

## SEC.21

2026-10-01 to 2026-10-01; 1 commits; not open.

- 2026-10-01: Lock out an address after failed logins (SEC.1, SEC.21) (#220)

## SEC.22

2026-10-01 to 2026-10-01; 1 commits; not open.

- 2026-10-01: Trusted devices, password blocklist and hashing, two-factor sign-in, fail2ban (SEC.22 to SEC.27) (#221)

## SEC.27

2026-10-01 to 2026-10-01; 3 commits; not open.

- 2026-10-01: TODO: bugs-first, then security plan; login brute-force research (SEC.20-SEC.27)
- 2026-10-01: TODO: always-on log (SEC.28) and fail2ban filter and jail examples (SEC.27)
- 2026-10-01: Trusted devices, password blocklist and hashing, two-factor sign-in, fail2ban (SEC.22 to SEC.27) (#221)

## SEC.28

2026-10-01 to 2026-10-01; 3 commits; not open.

- 2026-10-01: TODO: always-on log (SEC.28) and fail2ban filter and jail examples (SEC.27)
- 2026-10-01: TODO: sample log rotation for SEC.28 (logrotate, newsyslog, app-side on Windows)
- 2026-10-01: Always-on activity log and failed-login records (SEC.28, SEC.20) (#217)

## SEC.29

2026-10-08 to 2026-10-08; 2 commits; not open.

- 2026-10-08: Two-step sign-in on pyotp, QR codes on segno (SEC.29) (#612)
- 2026-10-08: Plan: retire SEC.29 and TEST.72; file TEST.103 and TEST.104 (#613)

## SEC.30

2026-10-08 to 2026-10-08; 2 commits; not open.

- 2026-10-08: Rate limits and login lockouts count on Redis; the lockout tests stop depending on the clock (SEC.30, TEST.83) (#643)
- 2026-10-08: Retire NAV.32, TEST.104, SEC.30, TEST.83; file TEST.105; fold GEN.73 into GEN.71; GEN.72 note (#644)

## SEC.31

2026-10-07 to 2026-10-07; 1 commits; not open.

- 2026-10-07: Bug fixes, lane 2 batch 1 (OPS.6, UX.36, MAP.117, NAV.41, OPS.7, SEC.31, UX.44, MAP.118) (#457)

## SEC.32

2026-10-23 to 2026-10-23; 2 commits; not open.

- 2026-10-23: SEC.32: Argon2id password hashing, a 1,024-character password limit and the __Host- session cookie (#2667)
- 2026-10-23: Docs: retire SEC.32 (#2667) (#2668)

## SEC.33

2026-10-21 to 2026-11-06; 2 commits; open.

- 2026-10-21: TODO: SEC.33 and TEST.261, the standing security and bug hunt audit release gates (Boss, 5:10 PM EDT) (#2342)
- 2026-11-06: Docs: file SEC.37 to SEC.42 from the 8.1 security first pass, log it on SEC.33, record TEST.658 (train tooling) (#5380)

## SEC.34

2026-10-24 to 2026-10-24; 3 commits; not open.

- 2026-10-24: docs: retire GEN.762, GEN.756, GEN.757, MAP.287, MAP.291; file SEC.34, SEC.35 and lane findings (#3133)
- 2026-10-24: SEC.34, SEC.35: Pillow 13.0.0 in both locks; the picture upload reads only PNG, JPEG and ICO
- 2026-10-24: docs: retire SEC.34, SEC.35; drop TEST.335 (#3141) (#3143)

## SEC.35

2026-10-24 to 2026-10-24; 3 commits; not open.

- 2026-10-24: docs: retire GEN.762, GEN.756, GEN.757, MAP.287, MAP.291; file SEC.34, SEC.35 and lane findings (#3133)
- 2026-10-24: SEC.34, SEC.35: Pillow 13.0.0 in both locks; the picture upload reads only PNG, JPEG and ICO
- 2026-10-24: docs: retire SEC.34, SEC.35; drop TEST.335 (#3141) (#3143)

## SEC.36

2026-10-24 to 2026-10-24; 3 commits; not open.

- 2026-10-24: docs: file SEC.36 (GIF and animated GIF avatars) (#3144)
- 2026-10-25: Account pictures: PNG, JPEG, ICO and GIF, animated too, checked under caps and stored as a clean copy (SEC.36)
- 2026-10-24: Docs: retire SEC.36, OPS.29, TEST.376, TEST.372, OPS.70; file TEST.379 (#3192)

## SEC.37

2026-11-06 to 2026-11-06; 4 commits; not open.

- 2026-11-06: Docs: file SEC.37 to SEC.42 from the 8.1 security first pass, log it on SEC.33, record TEST.658 (train tooling) (#5380)
- 2026-11-06: SEC.37: Boss chose option 1 (root-owned code tree) (#5414)
- 2026-11-06: SEC.37: the code tree is root-owned and read-only to the web user; config.json lives in the state folder (#5420)
- 2026-11-06: Retire SEC.37; file SEC.44 and TEST.661 (#5421)

## SEC.38

2026-11-06 to 2026-11-06; 2 commits; not open.

- 2026-11-06: SEC.38: the guard and the stop reasons do not trust the shared run folders (forged marker, planted link) (#5422)
- 2026-11-06: Retire SEC.38 after #5422 (#5425)

## SEC.39

2026-11-06 to 2026-11-06; 1 commits; not open.

- 2026-11-06: SEC.39: every queue carries JSON; a Redis user with a password is shipped and update.sh warns while Redis is open (#5437)

## SEC.40

2026-11-06 to 2026-11-06; 2 commits; not open.

- 2026-11-06: SEC.40, SEC.41: nltk 3.10.4 in both locks; the API refuses nan and inf in its numbers (#5451)
- 2026-11-06: Retire SEC.40 and SEC.41 (PR 5451) (#5452)

## SEC.41

2026-11-06 to 2026-11-06; 2 commits; not open.

- 2026-11-06: SEC.40, SEC.41: nltk 3.10.4 in both locks; the API refuses nan and inf in its numbers (#5451)
- 2026-11-06: Retire SEC.40 and SEC.41 (PR 5451) (#5452)

## SEC.42

2026-11-06 to 2026-11-06; 4 commits; not open.

- 2026-11-06: Docs: file SEC.37 to SEC.42 from the 8.1 security first pass, log it on SEC.33, record TEST.658 (train tooling) (#5380)
- 2026-11-06: File TEST.664; retire SEC.44; audit thread owns SEC.42 and SEC.45 (#5442)
- 2026-11-06: SEC.42: low-severity audit points (workflows, gunicorn unit, SVG check, hygiene) (#5455)
- 2026-11-06: Retire SEC.42 (PR 5455); OPS.110 scope per Boss (#5457)

## SEC.43

2026-11-06 to 2026-11-06; 4 commits; not open.

- 2026-11-06: File SEC.43: Update button in the admin pages (#5410)
- 2026-11-06: SEC.43 owner: the audit thread (#5413)
- 2026-11-06: SEC.43: Update button in the admin pages (#5424)
- 2026-11-06: Retire SEC.43 and GEN.1194 (#5426)

## SEC.44

2026-11-06 to 2026-11-06; 4 commits; not open.

- 2026-11-06: Retire SEC.37; file SEC.44 and TEST.661 (#5421)
- 2026-11-06: SEC.44: root only uses folders from config.json that pass a check (#5440)
- 2026-11-06: File TEST.664; retire SEC.44; audit thread owns SEC.42 and SEC.45 (#5442)
- 2026-11-06: SEC.44 follow-up: accept home folders in the config.json folder check (#5444)

## SEC.45

2026-11-06 to 2026-11-06; 3 commits; not open.

- 2026-11-06: File TEST.664; retire SEC.44; audit thread owns SEC.42 and SEC.45 (#5442)
- 2026-11-06: SEC.45: a worker runs only the functions the program queues (#5467)
- 2026-11-06: Retire SEC.45 (PR 5467); file TEST.667 (#5468)

## SEC.46

2026-11-06 to 2026-11-06; 1 commits; not open.

- 2026-11-06: Record SEC.46; retire TEST.663 and PERF.241 (#5447)

