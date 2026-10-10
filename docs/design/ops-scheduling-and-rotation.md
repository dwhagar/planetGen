# Update-script operations: Apache reload, daily maintenance, locks and backup rotation

How `update.sh` reloads the web server, how the daily
maintenance run is scheduled on Linux and macOS, which lock keeps
runs from overlapping, and how the settings JSON backups are rotated into 18
slots. The
reproducible-galaxy design these items belong to is
[reproducible-galaxies.md](reproducible-galaxies.md) (sections 6 to 8);
this note holds the measured behaviour, the commands and the tested code
the implementer needs.

Informs: OPS.8, OPS.13, OPS.14, OPS.16, OPS.17, OPS.18, ADM.19, ADM.20

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [S] seen in a search result or a fetched primary file (URL
in Sources), [C] computed or measured in the research sandbox, [R]
recalled and unconfirmed. The research environment could only read
search-result text and raw files, not whole papers or vendor manuals, and
had no macOS; every [R] below is on the
Evidence notes list for checking when access is allowed.

## Decisions already taken

- OPS.8 (Boss, 2026-10-01): "it should just automatically reload apache2
  if it's running as root." macOS (gunicorn) stays as it is
  unless Boss asks.
- OPS.16/OPS.18 (Boss, 2026-10-02 02:28Z): "a JSON file is only changed
  with the deltas at the end of the day. We're going to have to build a
  maintenance script for powershell and bash that will run the positional
  update script, then kick off this delta script that will update the JSON
  so it's only updated once every 24 hours, old JSON files are kept in the
  following order, 1 year ago, 6 months ago, 4 weeks ago, 7 days ago. A
  total of 18 backup slots."
- ADM.20 (Boss, 2026-10-02 02:31Z): merge now is Phase 3, low priority.

The uploaded source documents in this folder contain nothing on
scheduling, locking or rotation. The one overlap is `Library Migration
Workflow.md` and `Web UX Development Notes.md` (section 5), which list
APScheduler beside RQ; the daily run must survive a restart of the web
process, so the operating system's scheduler is used and RQ stays the
Generate queue.

## Summary

1. OPS.8 is simpler than written: the update script already requires
   root, so the "print the command" branch never runs. The
   closing message becomes a tested `reload_apache` step (section 1).
2. A graceful reload kills an in-flight mod_wsgi daemon request after about
   4 s; touching the WSGI file does not.
3. OPS.18 uses borg prune semantics on UTC names: exactly 18 distinct files,
   reaching 6 to 18 months back (restic's union semantics reach 4 months).
4. Locks are OS-level lock files (`filelock`), created and used by the web
   user, never PID files.
5. OPS.17 extends the monthly installers that already exist.

## 1. Reloading the web server from the update script (OPS.8)

### What the repo does today

`update.sh` ends by printing `sudo systemctl reload apache2` (`restart`
when `APACHE_NEEDS_RESTART` was set by `ensure_apache_modules` in
`scripts/deploy-common.sh`), or on macOS
`launchctl kill SIGHUP system/org.planetgen.gunicorn`. Behind nginx or Caddy on Linux the app
server is `planetgen-gunicorn`, which OPS.8 does not mention;
`examples/systemd/planetgen-update.service.d/reload-web.conf` already
reloads it from the timer. One function should cover both.

### Commands per platform

| Platform | Unit | Test config | Graceful | Hard restart |
|---|---|---|---|---|
| Debian/Ubuntu | `apache2` | `apache2ctl configtest` | `systemctl reload apache2` (`ExecReload=apachectl graceful`) | `systemctl restart apache2` |
| RHEL/Fedora/Rocky | `httpd` | `httpd -t` | `systemctl reload httpd` (`ExecReload=httpd $OPTIONS -k graceful`) | `systemctl restart httpd` |
| macOS Homebrew | brew service | `httpd -t` | `apachectl graceful` (unverified) | `brew services restart httpd` |

The Debian unit text was read from the installed Ubuntu 24.04 package [C];
the RHEL form is from third-party `httpd.service` mirrors [S]. The macOS
row was not run and matters little: macOS is nginx plus gunicorn.

### Graceful, restart and mod_wsgi daemon mode

Measured on Ubuntu 24.04, Apache 2.4.58, mod_wsgi 5.0.0, daemon mode
(`processes=1 threads=2`), with a 6 s request in flight [C]:

| Action | New daemon? | Config re-read? | In-flight request |
|---|---|---|---|
| edit a module file, no reload | no | no | n/a |
| `apache2ctl graceful` | yes | yes (a newly enabled `headers` module took effect) | 500 after about 4 s |
| same with `graceful-timeout=15` and `shutdown-timeout=15` | yes | yes | still 500 after about 4 s |
| `touch wsgi.py` (`WSGIScriptReloading On`, daemon mode) | yes, at next request | no | 200 after 6.0 s |
| `kill -USR1 <daemon pid>` with `graceful-timeout=15` | yes | no | 200 after 6.0 s |

Touching the script is the lightest reload and drains in-flight requests;
`graceful` re-reads Apache's config but does not drain mod_wsgi daemon
requests. This matches mod_wsgi's documentation: touching the script works in daemon
mode only, and `graceful-timeout` belongs to the daemon's own restart
signals, not to Apache's graceful restart [S]. planetGen's requests are
short (`request-timeout=60`), so reload is acceptable, but the docs should
say it can abort a request that is mid-flight and that `touch
src/html/wsgi.py` does not. A newly enabled module is loaded by a graceful
reload, so the "restart after `a2enmod`" rule in `update.sh` is stricter
than needed [C].

Generate jobs: `src/planetgen/web/job_runner.py` starts steps with
`start_new_session=True`, which leaves the daemon's session but not the
systemd cgroup. `apache2.service` has `KillMode=mixed`, so
`systemctl restart apache2` can kill a running Generate job and
`reload` does not [R; needs a systemd host to test].
[deployment/macos.md](../deployment/macos.md) and
`planetgen-gunicorn.service` already say "HUP leaves a running Generate
job alone". That is the strongest reason to prefer reload.

Recommended rule: reload always. Restart only if the script itself just
installed `libapache2-mod-wsgi-py3` and `jobs.active_job()` shows no
running Generate job; otherwise reload and print a note. The package's
postinst already tries `invoke-rc.d apache2 restart` [C], so a reload
afterwards is usually enough.

### Tested bash function

Run in four situations [C]: Apache stopped (rc 3, not started), running
with a good config (rc 0, daemon PID changed), running with a broken config
(rc 4, daemon PID unchanged, old server untouched), and the systemd path
with a stub `systemctl` (reload and restart verbs, rc 0). It is Bash 3.2
safe. Under `set -e` call it as `reload_apache || rc=$?`.

```bash
# Return codes: 0 reloaded; 3 not installed/not running (nothing done); 4 configtest failed
# (nothing done); 5 reload/restart failed or Apache is down afterwards.
_apache_unit() {                       # systemd unit running Apache, if any
    [[ -d /run/systemd/system ]] && command -v systemctl >/dev/null 2>&1 || return 1
    local u; for u in apache2 httpd; do
        systemctl cat "$u.service" >/dev/null 2>&1 && { echo "$u"; return 0; }
    done; return 1
}
_apache_ctl() { local c; for c in apache2ctl apachectl httpd; do command -v "$c" && return 0; done; return 1; }
_apache_running() {
    if [[ -n "$1" ]]; then systemctl is-active --quiet "$1.service"
    else pgrep -x apache2 >/dev/null 2>&1 || pgrep -x httpd >/dev/null 2>&1; fi
}
reload_apache() {
    local mode="${1:-reload}" unit ctl
    unit="$(_apache_unit)" || unit=""; ctl="$(_apache_ctl)" || ctl=""
    [[ -n "$unit$ctl" ]] || { echo "Apache: not installed."; return 3; }
    _apache_running "$unit" || { echo "Apache: not running; not starting it."; return 3; }
    if [[ -n "$ctl" ]] && ! "$ctl" configtest >/dev/null 2>/tmp/apache-configtest.$$; then
        echo "error: Apache's configuration test failed; not reloading:" >&2
        sed 's/^/  /' /tmp/apache-configtest.$$ >&2; rm -f /tmp/apache-configtest.$$; return 4
    fi
    rm -f /tmp/apache-configtest.$$
    if [[ -n "$unit" ]]; then systemctl "$mode" "$unit.service" || return 5
    elif [[ "$mode" == restart ]]; then "$ctl" restart || return 5
    else "$ctl" graceful || return 5; fi
    sleep 1; _apache_running "$unit" || { echo "error: Apache is down after the $mode." >&2; return 5; }
    echo "Apache: ${mode}ed (${unit:-$ctl})."
}
```

Notes. `apache2ctl graceful` on a stopped Apache starts it [C], hence the
running check first. The explicit configtest gives one clear message on
Debian and RHEL alike. A config that passes `-t` can still fail at run time
(port in use), hence the check that Apache is up afterwards. A real run
should also probe `curl -fsS --max-time 10 http://127.0.0.1/api/health`
through the vhost (the route [server-checklist.md](../server-checklist.md)
uses); this was not scripted because it needs the site config. The pull and
migration have already happened when the reload runs, so a failure cannot
roll back: print the manual command, set `RELOAD_FAILED=1` and end non-zero
so a scheduled unit shows as failed.

## 2. Scheduling the daily maintenance run (OPS.16, OPS.17)

### What already exists

`examples/maintenance/` has `install-maintenance-timer.sh` (systemd and
launchd),
`planetgen-update.{service,timer}`, `planetgen-orbits@.{service,timer}`
(one per database) and `examples/macos/org.planetgen.{orbits,update}.plist`.
All run monthly on the 1st (update 03:00, orbits 03:30) as root
with catch-up. The installer header says `install.sh` and `update.sh`
deliberately do not call it. So OPS.17 means: call the existing installers
from install and update, switch to daily, and retire the monthly orbit
timers (the daily run includes `planetgen.cli.orbits`; leaving both would
run orbits twice on the 1st). The header of `planetgen.cli.orbits` and
[architecture.md](architecture.md) call it "meant for a monthly timer";
both change to daily. The per-database instances mean OPS.16 and OPS.13
must iterate every galaxy database. A monthly `planetgen-update.timer` that
stays takes the same maintenance lock and waits up to 30 minutes for it.

### Defaults per OS

| | Default | Catch-up after downtime |
|---|---|---|
| Linux with systemd | systemd timer plus `Type=oneshot` service | `Persistent=true` runs once at next start if a trigger was missed |
| Linux without systemd (containers, Alpine) | `/etc/cron.d/planetgen-maintenance` | none; add `@reboot` and `--if-due` |
| macOS | launchd daemon, `StartCalendarInterval` | asleep: runs at wake; powered off: skipped; add `RunAtLoad` and `--if-due` |

anacron is not the default (whole days only, root only, usually absent on
servers [S]). In `/etc/cron.d` the file must be root-owned, not group or
other writable, carry a user field, and have no dot in its name or Debian
cron skips it silently [S]. Detect systemd with `[ -d /run/systemd/system ]`
and `systemctl` on the path.

Run as the web user on Linux and macOS, not root: the web app writes the
settings JSON and Generate job state, root-run maintenance would leave
root-owned files the web user cannot rotate, and a root-created lock file
locks the web user out (section 3). The code tree stays root-owned
(`set-permissions.sh`), so a web-user process cannot change code that root
later runs. `config.json` is `root:<apache group>` mode 640, readable by the
web user, and `detect_apache_group` in `examples/apache/apache-identity.sh`
finds the account. The monthly orbit timer runs as root with
`/etc/planetgen/maintenance.env`; that is what is retired.

Implement the work once in Python as `planetgen.cli.maintenance` (steps,
lock, `--if-due`, `--dry-run`, `--only STEP`, rotating log, exit codes
below) and keep `scripts/maintenance.sh` as a thin launcher that
finds the right Python and runs it.

### Linux: systemd timer

`systemd-analyze verify` on systemd 255 accepted both units without
warnings and flagged a deliberately misspelled key [C]. `systemd-analyze
calendar '*-*-* 03:30:00 UTC'` printed the expected elapses [C].

```ini
# /etc/systemd/system/planetgen-maintenance.service
[Unit]
Description=planetGen daily maintenance (orbits, settings merge, backup rotation)
After=network-online.target mysql.service mariadb.service redis-server.service
Wants=network-online.target
StartLimitIntervalSec=6h
StartLimitBurst=3

[Service]
Type=oneshot
User=www-data
Group=www-data
WorkingDirectory=/var/lib/planetGen
ExecStart=/var/lib/planetGen/scripts/maintenance.sh --if-due
SuccessExitStatus=75          # EX_TEMPFAIL: lock held or Generate job running; not a failure
Restart=on-failure
RestartSec=15min
TimeoutStartSec=6h
Nice=10
IOSchedulingClass=idle
NoNewPrivileges=true

# /etc/systemd/system/planetgen-maintenance.timer
[Unit]
Description=Run planetGen daily maintenance
[Timer]
OnCalendar=*-*-* 03:30:00 UTC
Persistent=true
AccuracySec=1min
[Install]
WantedBy=timers.target
```

- `Persistent=true` stores the last trigger time and fires at once if one
  was missed while the timer was inactive; it works only with `OnCalendar=`
  [S]. `systemctl clean --what=state planetgen-maintenance.timer` removes
  the stamp on uninstall [S]; the stamp path
  `/var/lib/systemd/timers/stamp-<unit>.timer` is [R].
- No `RandomizedDelaySec`: the existing timers use fixed times so ordering
  is guaranteed.
- UTC in `OnCalendar` keeps the run's UTC date stable (backup names are
  UTC); a local schedule can put two runs in one UTC date or none across a
  DST change, which rotation tolerates. launchd offers
  local time only.
- Test: `sudo systemctl start planetgen-maintenance.service`,
  `systemctl list-timers 'planetgen-*'`, `journalctl -u planetgen-maintenance`.
  Output goes to journald.

Idempotent create, leave, refresh, remove:

- Create: install the two units with `install -m 644` only if the `.timer`
  file does not exist; `systemctl daemon-reload`; `systemctl enable --now
  planetgen-maintenance.timer` (a no-op when already enabled).
- Leave: if the timer file exists, do nothing, including after an admin ran
  `systemctl disable` (that is how they turn it off). The schedule changes
  by drop-in (`systemctl edit`, with an empty `OnCalendar=` line first to
  clear the default); update never touches drop-ins.
- Refresh: ship a `# planetgen-unit-version: N` line and rewrite only the
  `.service` when N differs; the `.timer` (the schedule) is left alone.
- Remove: `systemctl disable --now planetgen-maintenance.timer`, delete both
  files, `daemon-reload`.
- A deleted file cannot be told from a never-installed one, so give an
  off switch that survives updates: `maintenance.schedule` in `config.json`
  (`"auto"` default, `"off"`), read by install and update through
  `planetgen.util.appconfig`. It works on Linux and macOS.

### Linux without systemd: cron

```
# /etc/cron.d/planetgen-maintenance   (root:root 0644; no dot in the name)
SHELL=/bin/bash
PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
30 3 * * * www-data /var/lib/planetGen/scripts/maintenance.sh --if-due >> /var/log/planetgen/maintenance.log 2>&1
@reboot    www-data sleep 300 && /var/lib/planetGen/scripts/maintenance.sh --if-due >> /var/log/planetgen/maintenance.log 2>&1
```

Cron uses local time and has no catch-up; `--if-due` makes the `@reboot`
line safe.

### macOS: launchd

Apple's scheduling guide says a `StartCalendarInterval` job whose time
passes while the Mac is asleep runs at wake; if the Mac is off it waits for
the next scheduled time [S]. Coalescing of several missed intervals is
documented for `StartInterval`; for `StartCalendarInterval` no statement
was found. [deployment/macos.md](../deployment/macos.md) says only "a run
missed while the Mac was switched off is skipped"; it should add that an
asleep Mac runs at wake. `RunAtLoad` plus `--if-due` covers the powered-off
case; `pmset repeat wakeorpoweron` can schedule a wake [R].

```xml
<!-- /Library/LaunchDaemons/org.planetgen.maintenance.plist (parses with plistlib; xmllint-clean) -->
<plist version="1.0"><dict>
  <key>Label</key><string>org.planetgen.maintenance</string>
  <key>ProgramArguments</key><array>
    <string>/bin/bash</string><string>/var/lib/planetGen/scripts/maintenance.sh</string><string>--if-due</string></array>
  <key>UserName</key><string>_www</string>
  <key>WorkingDirectory</key><string>/var/lib/planetGen</string>
  <key>StartCalendarInterval</key><dict><key>Hour</key><integer>3</integer><key>Minute</key><integer>30</integer></dict>
  <key>RunAtLoad</key><true/>
  <key>LowPriorityIO</key><true/><key>Nice</key><integer>10</integer>
  <key>StandardOutPath</key><string>/usr/local/planetgen/log/maintenance.log</string>
  <key>StandardErrorPath</key><string>/usr/local/planetgen/log/maintenance.log</string>
</dict></plist>
```

`install_launchd_plist` in `scripts/deploy-common.sh` already leaves an
existing plist alone and otherwise copies, chowns root:wheel, sets 644 and
runs `launchctl bootstrap system ... || launchctl load -w`; reuse it
(`bootstrap` fails if the label is loaded, so the "present, leave it"
check must come first, as the helper does). Remove with `sudo launchctl
bootout system/org.planetgen.maintenance` and delete the plist. Test with
`sudo launchctl kickstart -k system/org.planetgen.maintenance`, `plutil
-lint`, `launchctl print`. Add a `newsyslog` rule for `maintenance.log`
beside the existing debug-log rule. launchd could not be run here; the
plist parses [C] but is untested on a Mac.

### `--if-due` and exit codes

Schedulers differ (Persistent catch-up, `RunAtLoad`, manual starts, merge
now), so the script decides whether a run is due: it records the last
successful completion (UTC) in the control database or a state file beside
the lock, and with `--if-due` exits 0 quickly if that was under about 20
hours ago. This enforces Boss's "only once every 24 hours" however often
the scheduler fires. Without `--if-due` it always runs (manual runs and
tests).

| Code | Meaning |
|---|---|
| 0 | done, or nothing due |
| 1 | a step failed |
| 2 | usage error |
| 75 | skipped: another maintenance or update run holds the lock, or a Generate job is running (EX_TEMPFAIL) |

Systemd gets `SuccessExitStatus=75`; for launchd the
code is informational. Alert when the run was skipped N days in a row, so a
multi-day Generate run cannot starve maintenance silently.

## 3. Lock files (OPS.16, ADM.20)

### What exists

`src/planetgen/web/jobs.py` holds the one-running-Generate-job lock by
hand: an `active` file created atomically (temp file plus `os.link`,
fallback `O_EXCL`), stale handling that looks up the named job and checks
the runner's PID (`/proc/<pid>/cmdline` on Linux to catch PID reuse), and a `.clearing` guard
so two admins cannot both clear a stale lock. That shape suits a lock that
must carry an owner id for the UI. The repo has no `filelock`,
`portalocker` or `fcntl` use outside a test helper
(`src/tests/worker_patches.py`). `redis` and `croniter` are already
dependencies.

### Options

| | Stale locks | Python 3.9 | Verdict |
|---|---|---|---|
| `fcntl.flock` by hand | none | n/a | a few lines, but the corner cases are yours |
| `filelock` 4.0.12 (MIT) | none for the hard locks | needs 3.10+; last 3.9 release 3.19.1 | **recommended**, pin by marker |
| `portalocker` 4.4.0, `fasteners` 0.20 | none | 3.2.0 / any | alternatives |
| PID file only | needs a PID check, races with PID reuse | n/a | avoid |
| Redis key with TTL | TTL | n/a | avoid: Redis may be down exactly when maintenance must run |

Versions read from the PyPI JSON API [C]. `requirements.lock` already uses
per-Python markers, so `filelock==3.19.1 ; python_full_version < '3.10'`
and the current release for 3.10+ fit `scripts/lock-requirements.sh`.

### Measured with `filelock` 3.19.1 on Linux [C]

1. A second acquirer with `timeout=0` raises `Timeout` at once (0.000 s).
2. `kill -9` of the holder frees the lock immediately; the file stays on
   disk, harmlessly.
3. A lock file created by root with mode 0644 fails for user `nobody` with
   `PermissionError: [Errno 13]`, not `Timeout`. `mode=0o666` (or 0o660 and
   a shared group) works. Create and use the lock as one user and have the
   web code catch `PermissionError`.
4. `/proc/sys/fs/protected_regular` was 0 here; at 1 or 2, opening a file
   you do not own with `O_CREAT` in a sticky world-writable directory such
   as `/tmp` fails even for root [R]. Keep locks in a `locks/` directory
   under the data directory (path from `config.json`).

Never store the holder's PID inside the lock file;
put diagnostics in a sidecar written atomically (`maintenance.lock.json`:
pid, started UTC, host, command). Two lock objects in one process
conflict (`flock` belongs to the open file description), which is good for a
second merge-now click. Locking over NFS or SMB is unreliable [R]; keep the
data directory on a local disk.

### What to lock

| Lock | Held by | Purpose | If it cannot be taken |
|---|---|---|---|
| `locks/maintenance.lock` | `planetgen.cli.maintenance`, `update.sh`, merge now (ADM.20) | no two of: daily run, update, merge now | maintenance: exit 75 and log "skipped"; update: wait up to 30 min, then stop with a message; merge now: show "the daily run holds the lock" and do nothing |
| `jobs/active` (existing) | Generate page and API | one Generate job | maintenance checks `jobs.active_job()` and exits 75; `start_job` refuses while the maintenance lock is held (try-lock, release) |

No concrete corruption case was found for maintenance beside a Generate
job (orbit update rewrites positions in lockstep; the merge clears pending
deltas), but "skip and run tomorrow" is free: the orbit update catches up
to the current time and the merge carries more deltas. The delta merge
(GEN.61) should clear pending rows with `id <= highest id read`, never
"everything", because admin edits do not take the lock.

Order in the update script: take the lock first, before the `git reset
--hard` that changes code a running maintenance process imports. `update.sh`
runs on the system Python before step 2 checks libraries, so on the first
update that adds `filelock` the import fails: either warn and carry on
without the lock for that one run, or run a stdlib-only helper
(`"$PYTHON" scripts/runlocked.py ...`). macOS has no `flock(1)` [R].

## 4. Backup rotation: 18 slots, 7 daily, 4 weekly, 6 monthly, 1 yearly (OPS.18)

### What gets rotated

Files named `<32 hex seed>-<22 hex key>-<YYYYMMDD>-<HHMMSS>Z.json` in the
data directory (ADM.18, [reproducible-galaxies.md](reproducible-galaxies.md)
section 7). A new file is written only when deltas are pending or the epoch
changed, so quiet days leave gaps, and merge now (ADM.20) can add several
files in one day. Rotation must not assume one file per day.

### Existing implementations

borg prune [S] walks each tier (day, ISO week, month, year) newest to
oldest, keeps the newest file of each new period, and does not count a file
a finer tier already kept. restic forget [S] ORs its rules, so the newest
snapshot satisfies all four tiers. rotate-backups [S] keeps the oldest of
each slot by default. Borg 1.2+ also keeps the oldest archive if a rule fell
short; that would add a 19th file, so leave it out.

### Comparison by simulation [C]

Three years of simulated merges, rotation after every day, counts 7/4/6/1:

| Scenario | Rule | Files kept | Oldest file |
|---|---|---|---|
| a file every day | borg-style | 18, never more | 369 days |
| a file on 30% of days | borg-style | 18 | 368 days |
| 1 to 3 files per day | borg-style | 18 | 368 days |
| a file every day | restic-style (union) | 13 (max 15 during the run) | 126 days |
| a file on 30% of days | restic-style | 12 | 127 days |
| 1 to 3 files per day | restic-style | 13 | 125 days |

Restic-style fails Boss's goal ("1 year ago, 6 months ago, 4 weeks ago, 7
days ago"): its yearly slot is the newest file of the current year, the same
file as the daily slot, so history reaches only about four months.
Borg-style makes the tiers hold distinct files and walks each into older
periods. Over 730 consecutive daily files it kept exactly 18 on every day;
the newest weekly file was 7 to 13 days old, the oldest weekly 28 to 34, the
oldest monthly 157 to 218, and the yearly file 188 to 553 days old (always
the last file of a past calendar year). That last range is the honest "1
year ago": a rolling "age 365" target would need every file kept until its
year is up, because a file deleted at age 15 days is gone by age 365.
Calendar periods are stable under repeated rotation; rolling targets are
not.

Result for 2026-10-09 with a file every day since 2024-01-01 and rotation
after each:

| Tier | Kept (age in days) |
|---|---|
| daily (7) | 2026-10-09 (0), 10-08, 10-07, 10-06, 10-05, 10-04, 10-03 (6) |
| weekly (4) | 09-27 (12), 09-20 (19), 09-13 (26), 09-06 (33) |
| monthly (6) | 09-30 (9), 08-31 (39), 07-31 (70), 06-30 (101), 05-31 (131), 04-30 (162) |
| yearly (1) | 2025-12-31 (282) |

### Algorithm

1. Parse names matching `^<32 hex>-<22 hex>-<YYYYMMDD>-<HHMMSS>Z\.json$`;
   ignore everything else (downloads, `.tmp` files, other galaxies). Group
   by seed: one rotation per galaxy.
2. The time is the UTC time in the name. Never mtime (copies, restores and
   antivirus change it) and never local time.
3. For each tier in (daily, weekly, monthly, yearly) with `n` = 7, 4, 6, 1:
   walk files newest to oldest; when a file's period (UTC date; ISO year
   and week, Monday start; UTC year-month; UTC year) differs from the
   previous file's, mark it seen and, if the file is not already kept, keep
   it and count it; stop at `n` counted.
4. The newest file is the first daily keeper, so the current file is always
   kept.
5. Several files in one day or week: the newest of the period wins. A merge
   now file is superseded by the night run, as ADM.20 expects.
6. Missed days cost nothing: only periods that contain a file count.
7. A file wanted by two tiers is kept once, labelled with the finest tier,
   and the coarser tier moves on to the next period. That is why the total
   is exactly 18 once history exists.
8. Fewer than 18 files: keep everything the rules select.

DST is irrelevant because periods use UTC names. Leap days and the ISO
week-year boundary (2025-12-29 is in 2026-W01; 2027-01-01 is in 2026-W53)
are covered by tests. Clock skew has two guards: the merge names its file
`max(now, newest_existing + 1 s)` so a stepped-back clock cannot produce a
name older than the current file, and rotation deletes nothing, and warns,
if any name is more than a day in the future (one future-dated file would
become "newest" and push the real backups into the weekly and monthly
slots).

Delete only after the new file is written, fsynced and read back, under the
maintenance lock, with a `--dry-run` that prints "keep (slot) / delete".

### The pure function

The function never reads the clock, so tests inject dates. It uses no
syntax newer than Python 3.8; tests ran on 3.11 and 3.13 only (no 3.9
interpreter was available).

```python
from collections import OrderedDict
from datetime import timedelta

TIERS = ("daily", "weekly", "monthly", "yearly")
DEFAULT_COUNTS = OrderedDict([("daily", 7), ("weekly", 4), ("monthly", 6), ("yearly", 1)])

def _period(tier, when):                      # `when` is an aware UTC datetime
    if tier == "daily":   return (when.year, when.month, when.day)
    if tier == "weekly":  iso = when.isocalendar(); return (iso[0], iso[1])   # ISO year, week
    if tier == "monthly": return (when.year, when.month)
    return (when.year,)

def select_keep(stamps, counts=DEFAULT_COUNTS):
    """{timestamp: tier} for the files to keep, newest first (borg prune semantics)."""
    ordered = sorted(set(stamps), reverse=True)
    kept = OrderedDict()
    for tier in TIERS:
        n, last, taken = counts[tier], None, 0
        if n <= 0:
            continue
        for ts in ordered:
            p = _period(tier, ts)
            if p == last:
                continue
            last = p                         # period seen, even if its newest file was kept by a finer tier
            if ts in kept:
                continue
            kept[ts] = tier
            taken += 1
            if taken == n:
                break
    return OrderedDict(sorted(kept.items(), reverse=True))

def next_stamp(existing, now):               # name for a new file: never at or before the newest
    now = now.replace(microsecond=0)
    return max(existing) + timedelta(seconds=1) if existing and now <= max(existing) else now

def suspicious_future(stamps, now, slack=timedelta(days=1)):
    return sorted(ts for ts in stamps if ts > now + slack)
```

The prototype also had `parse_name`, a `mode="union"` switch reproducing
restic for the comparison, and `to_delete`.

### Tests (15, all passing on 3.11 and 3.13) [C]

Name parsing (the ADM.18 example parses; `config.json` and an impossible
date do not); fewer files than slots deletes nothing; three years of one
file per day with daily rotation stays at most 18, ends exactly 18 split
7/4/6/1, keeps the newest, and the oldest kept is the last file of the
previous calendar year; rotating only every 9th day stays valid; rotating
a rotated set changes nothing; 300 random sets (30 days to 7 years) keep
the newest, at most 18, idempotent; several files in one day keep only the
last once history is long enough; a 40-day gap leaves the 7 daily slots on
the 7 newest days with a file; the ISO week-year boundary (2025-12-29 and
2026-01-04 share a week; 2027-01-01 is 2026-W53); leap day 2028-02-29; UTC
day, not local day (23:30Z and 00:30Z the next day differ); no missing or
doubled day across a DST change; `next_stamp` (same second goes +1 s, a
stepped-back clock goes after the newest); the future-file guard flags a
2030 file against a 2026 clock but not one within a day of skew, and a
future-dated file without the guard takes a daily slot.

Put them in `src/tests/test_backup_rotation.py` (no filesystem), plus one
temp-directory test that creates files, rotates, and checks that
unparseable names and other seeds survive.

## 5. OPS.13 history table and atomic writes

### History table gotchas (keep last 10)

Run on MariaDB 10.11.14 only; MariaDB 11.4 and MySQL 8.4 were not tested [C]:

- `DELETE ... WHERE id NOT IN (SELECT id ... ORDER BY id DESC LIMIT 10)`
  fails ("doesn't yet support LIMIT & IN/ALL/ANY/SOME subquery"; MySQL has
  the same restriction [R]). Wrap it in a derived table: `NOT IN (SELECT id
  FROM (SELECT id FROM version_history WHERE galaxy_id=1 ORDER BY id DESC
  LIMIT 10) k)`. A cutoff form also works: `id < (SELECT id FROM (SELECT id
  ... ORDER BY id DESC LIMIT 1 OFFSET 9) k)`.
- Order and trim by the auto-increment `id`, not by date, so a wrong clock
  cannot evict newer rows. Store the date as an explicit UTC `DATETIME(6)`,
  not a `TIMESTAMP`.
- Insert and trim in one transaction, after the update succeeded, under the
  maintenance lock.
- A repeated update with no change would fill the 10 slots with identical
  rows and push out the older distinct versions ("Already up to date" is
  common). Record a row only when the key or any stored hash differs from
  the newest row; otherwise touch a `last_seen` column. The
  TODO's test (11 runs, oldest row goes) must use 11 distinct tuples, with a
  second test that 11 identical runs keep one row.
- "Per galaxy" means every galaxy database (`PLANETGEN_MYSQL_DATABASE_PREFIX`
  installs exist); the step iterates the galaxies in the control database.
- Line endings change hashes. `.gitattributes` forces LF only for `*.sh` and
  `src/html/**/*.py`. `requirements.lock`, `requirements-server.lock` and
  `src/planetgen/names/offensive_words.txt` are tracked text; a
  checkout with `core.autocrlf=true` gets CRLF and a different SHA-256 for
  identical content, so OPS.14 would warn "changed" between two servers
  for no real reason [C]. Fix with `*.lock text eol=lf` and
  `src/planetgen/names/*.txt text eol=lf`, or hash normalised text (read as
  text, normalise `\r\n`, encode UTF-8) in sorted file order.

### Atomic file replacement

- POSIX: temp file in the same directory, flush, `fsync`, `os.replace`,
  then `fsync` the directory. A hot reader thread saw 201,721 reads and 0
  partial reads over 2,000 replaces [C].
- macOS: `os.fsync` does not flush the drive cache; `fcntl.F_FULLFSYNC`
  does [R].
- The settings JSONs are never replaced: each merge writes a new name, so
  the primitive is "create atomically, never clobber". `os.link(tmp, final)`
  then `os.unlink(tmp)` raises `FileExistsError` if the name exists (as
  `jobs._write_lock` does); filesystems without hard links need an
  `O_EXCL` fallback. A second write to the same name was refused and left
  the content unchanged [C], and `next_stamp` gives same-second merges
  distinct names.
- Read the file back and compare before clearing pending rows. Clean
  leftover `.settings-*.tmp` files older than a day at the start of a run.

## Evidence notes

Measured [C] on Ubuntu 24.04, Apache 2.4.58, mod_wsgi 5.0.0, systemd 255
(not PID 1), MariaDB 10.11.14, Python 3.11 and 3.13: graceful versus touch
versus SIGUSR1; configtest; module load by graceful; the `reload_apache`
paths (systemd path with a stub `systemctl`); rotation simulations and
tests; `filelock` behaviour; the DELETE forms; the plist XML error; the
`os.replace` reader test; unit-file verification.

Still [R], to verify when paper and platform access is allowed: RHEL
`apachectl` and Red Hat's own `httpd.service`; macOS Apache graceful
forms; the systemd cgroup kill of
Generate jobs on `restart`; the systemd stamp path and
`StandardOutput=append:`; launchd coalescing and `pmset repeat`; `F_FULLFSYNC`; `protected_regular`; NFS and SMB `flock`; no `flock(1)` on
macOS; MySQL 8.4 matching MariaDB for the DELETE forms.

## Sources

- mod_wsgi: https://modwsgi.readthedocs.io/en/master/user-guides/reloading-source-code.html ; https://modwsgi.readthedocs.io/en/latest/release-notes/version-4.1.0.html
- Debian bug on `apache2ctl graceful`: https://bugs.edge.launchpad.net/bugs/1832182 ; `httpd.service` mirror: https://gitcode.com/src-openeuler/httpd/blob/master/httpd.service
- systemd.timer: https://manpages.debian.org/buster/systemd/systemd.timer.5 ; Debian cron(8): https://manpages.debian.org/jessie/cron/cron.8.en.html ; anacrontab(5): https://unix.com/man-page/linux/5/anacrontab/
- Apple, Scheduling Timed Jobs: https://developer.apple.com/library/mac/documentation/MacOSX/Conceptual/BPSystemStartup/Chapters/ScheduledJobs.html
- borg: https://borgbackup.readthedocs.io/en/1.4.1/usage/prune.html and `src/borg/archiver.py` on the 1.4-maint branch ; restic: https://restic.readthedocs.io/en/v0.12.0/060_forget.html ; rotate-backups: https://rotate-backups.readthedocs.io/en/latest/readme.html
- PyPI JSON API for filelock, portalocker, fasteners
