# Login brute-force protection

Research notes and the recommended design for protecting planetGen's
login against password guessing. Written 2026-10-01 for Boss's request
to plan the login blocking work before any code changes. Which pieces
are built, and in what order, is tracked in `docs/TODO.md`.

Status: every step is built (SEC.20, SEC.1, SEC.21, SEC.23, SEC.22, SEC.24,
SEC.25, SEC.26 and SEC.27; the activity log of SEC.28 carries the log
lines). Two-factor sign-in is optional per admin; a trusted device does
not skip the code; only the command line resets another admin's.
Passwords are now hashed with PBKDF2-SHA256 at 600,000 rounds (about
175 ms on the build machine, next to no memory), re-hashed on login.
Section 1 describes the site as it was before them.

## 1. What the site had before (2026-10-01)

- **Per-IP rate limit.** `POST /api/auth/login` allows 10 attempts a
  minute per client address (`auth.LOGIN_RATE_LIMIT`, Flask-Limiter,
  `memory://` storage). The web `/login` form calls the same route
  in-process and passes the real visitor's address through.
- **Per-username backoff.** `src/planetgen/api/loginbackoff.py`: 10 free
  failures per username, then a lock of 1 s doubling to 15 minutes,
  answered with 429 and `Retry-After` before the password is checked.
  Unknown usernames are counted like real ones, so a lock doesn't reveal
  which usernames exist. It is kept in memory: a restart forgets it, and
  each worker process keeps its own counts (the deployment guides run
  one process with five threads, so today there is one set of counts).
- **Passwords.** At least 12 characters, not equal to the username,
  hashed with werkzeug's `generate_password_hash` (werkzeug 3.1.9 from
  the lock files; its default is scrypt with N = 2^15, r = 8, p = 1).
  Login time is the same for known and unknown usernames.
- **Sessions.** `HttpOnly`, `Secure`, `SameSite=Strict` cookie, a fixed
  lifetime, stored as a SHA-256 hash; changing credentials ends every
  other session.

Gaps found while reading the code:

1. **No record of failed logins.** Nothing writes a failed or locked
   login anywhere (the audit log only records admin actions), so there
   is nothing to alert on and nothing for an external tool to watch.
2. **The web form hides failures from the access log.** A wrong
   password on `/login` answers 200 with an error message
   (`web/admin_pages.py`), so Apache's access log can't tell a failed
   login from a page view.
3. **The current-password check on `/account` isn't throttled.**
   `POST /api/auth/change-credentials` re-checks the current password,
   but a wrong one isn't counted by the backoff, and the web `/account`
   form only falls under the shared 300-a-minute page limit (in-process
   API calls skip the API's default limits). Someone holding a stolen
   session cookie can guess the password there quickly.
4. **Lockouts are easy to trigger on purpose.** Anyone who knows an
   admin's username can keep it locked for 15 minutes at a time. The
   planned per-IP lockout (up to a day) adds the same risk for everyone
   behind a shared address, or for every visitor at once if the site
   sits behind a proxy whose address isn't unwrapped (`proxy_fix`).

## 2. Standard methods

| Method | Stops | Costs | Fit here |
|---|---|---|---|
| Per-IP throttling | One machine guessing fast | Blocks everyone behind a shared address (NAT, office, proxy) | Yes, already built (10 a minute); the lockout adds a longer penalty |
| Per-account throttling with exponential backoff | Slow guessing of one account from many addresses | Lockout-as-DoS against the real owner | Yes, already built |
| Hard lockout (account disabled until an admin unlocks) | Everything | An attacker can disable every account | No: delays that expire are the standard choice |
| Trusted-device cookie | Lockout-as-DoS: a browser that has logged in before skips the per-account lock | A signed cookie and a small table | Yes, recommended |
| CAPTCHA after N failures | Bots | Third-party service, privacy, accessibility; solvers are cheap | No: needs an outside service on a self-hosted site |
| Proof-of-work challenge | Bulk guessing | JavaScript, CPU on phones | No: the throttles already bound guessing |
| Breached and common password blocklist | Credential stuffing and dictionary guesses | A word list in the repo, checked when a password is set | Yes, recommended |
| Slow password hash | Offline cracking after a database leak | CPU and memory per login | Yes, already in place; check the cost |
| Two-factor (TOTP) | A guessed or leaked password alone | Authenticator app, recovery codes | Yes, the strongest single step, for admins |
| Logging failed logins with the address | Nothing alone; enables alerts and outside tools | One log line and one table | Yes, first |
| fail2ban watching a log | Repeat offenders, at the firewall before Apache | Root setup on the server, Linux only | Yes, optional, documented |
| mod_evasive or mod_security at Apache | Floods and known attack patterns | Tuning, false positives | Mention only |

What the references say:

- **NIST SP 800-63B-4** (2025): the verifier must rate-limit failed
  attempts on an account, at most 100 consecutive failures, and may use
  delays that grow, CAPTCHA, or address allowlisting; it prefers delays
  over permanent lockout. It asks for a minimum of 15 characters for a
  password used alone (8 when it is one of two factors), support for at
  least 64 characters, no composition rules, and a check against a list
  of breached, common and expected passwords.
- **OWASP** (Authentication, Credential Stuffing and Password Storage
  cheat sheets): throttle both per account and per address; prefer
  lockouts that expire and grow over permanent ones; let a known device
  through a per-account lock ("device cookies"); keep error messages the
  same for every failure; log every failure with its address; hash with
  argon2id (19 MiB, 2 passes), scrypt (N = 2^17, r = 8, p = 1) or
  PBKDF2-HMAC-SHA256 (600,000 rounds); offer multi-factor sign-in.

## 3. Recommended design

In build order. Each step stands on its own.

1. **Log every failed and locked login.** One line per failure in the
   debug log (or a dedicated auth log) with the time, the client address
   and the username, and a row in the control database's audit log
   (`login.failed`, `login.locked`). The web `/login` form answers 401
   for a wrong password so the access log shows it too. This is what
   alerts, the admin page and fail2ban read.
2. **Per-IP lockout in Redis** (Boss's numbers; first in the control database's `login_throttle` table, moved to Redis by SEC.30): 3
   failures from one address lock it for 5 minutes, doubling with each
   further lockout up to 1 day; 429 and `Retry-After`, checked before the
   password, shared by every worker and kept across restarts. Safeguards
   that come with it:
   - `127.0.0.1`, `::1` and a configured allowlist are never locked.
   - If every request comes from one address that looks like a proxy
     (Apache in front with `proxy_fix` unset), the app warns at start-up
     and in the admin stats page instead of locking that address.
   - A successful login from an address resets its failure count; the
     doubling level decays (for example halves after a day without a
     lockout) instead of resetting at once, so an attacker can't clear it
     with one right guess.
   - An admin command (`planetgen` or a small `planetgen-admin`
     script) and an admin page list and lift lockouts.
3. **Move the per-username backoff into the same store**, so it is
   shared across workers and restarts; the numbers stay as they are.
4. **Trusted-device cookie.** A successful login sets a long-lived,
   signed, per-account device cookie (its hash stored in the control
   database). A login from a browser holding a valid device cookie for
   that username skips the per-username lock (the per-IP limit still
   applies), so an attacker can't lock the real admin out. Changing
   credentials revokes the account's device cookies.
5. **Count wrong current passwords.** A wrong current password on
   change-credentials counts as a failed login for that username and
   address, with the same lock.
6. **Password blocklist.** When a password is set, refuse it if it is on
   a bundled list of common and breached passwords (for example the top
   100,000 from a public breach corpus, shipped in the repo, offline).
   An online check against Have I Been Pwned's k-anonymity range API is
   an option that sends only the first five characters of the hash.
7. **Hashing cost.** Measure the current scrypt cost on the target
   server; if a login takes well under 250 ms, raise it (scrypt N = 2^16
   or 2^17, or PBKDF2-SHA256 at 600,000 rounds if memory is tight, since
   scrypt N = 2^17 uses 128 MiB per login). Re-hash on the next
   successful login when a stored hash uses older settings.
8. **Two-factor sign-in (TOTP) for admins.** Optional per account at
   first, with single-use recovery codes stored hashed; the Owner (once
   user accounts exist) may be required to use it. The second step is
   throttled like the password step.
9. **fail2ban recipe** in `docs/deployment/` (Linux): a filter for the
   auth log line from step 1 and a jail that bans an address at the
   firewall after repeated failures, with the server's own addresses
   ignored. Optional, not installed by `install.sh`.

Alerts (an email to admins after a lockout, or many failures) wait for
SMTP settings in the user-accounts work; until then the admin stats page
shows recent failures and current lockouts.

## 4. Decisions for Boss

Defaults are taken so work can start; each can change.

- Keep both limits (per address and per username)? Default: yes.
- Does a correct password from a locked address get in? Default: no,
  the lock is checked first, as the username backoff does now.
- Minimum password length: stay at 12, or 15 as NIST asks for a
  password used alone? Default: 12 plus the blocklist until two-factor
  exists, then revisit.
- Is the blocklist offline only, or may the server call Have I Been
  Pwned? Default: offline only.
- Two-factor optional or required for admins? Default: optional, and
  required for the Owner once roles exist.

## 5. Planned: libraries (2026-10-07)

The design above stays; its implementation moves to libraries
(library-migration.md): request limits to Flask-Limiter with Redis
storage, one-time codes to pyotp and QR codes to segno. The lockout
rules, the activity log and the fail2ban lines don't change.
