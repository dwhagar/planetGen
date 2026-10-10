# Blocking login brute force with fail2ban

planetGen locks an address out by itself after failed sign-ins (3
failures, then 5 minutes doubling to a day; [`api.md`](../api.md)). On
Linux, fail2ban can go further and drop a guessing address at the
firewall, so it can't reach the site at all. It reads the activity log
([`config.md`](../config.md#the-activity-log)), which writes one line per
failed sign-in with the client's address.

The filter and jail are in [`examples/fail2ban/`](../../examples/fail2ban/).
They match four events:

| Event | When |
|---|---|
| `AUTH login.failed` | A wrong username or password at `/login` or `POST /api/auth/login`. |
| `AUTH login.locked` | A sign-in refused, unchecked, while planetGen's own lockout is on. |
| `AUTH password.failed` | A wrong current password on the account page (change credentials). |
| `AUTH totp.failed` | A wrong two-factor code at sign-in, or when turning two-factor sign-in off. |

## Before you start

- The activity log must be written. `update.sh` creates
  `/var/log/planetgen/` and the log (`planetgen-log`); check that it
  fills when you fail a sign-in on purpose:
  `sudo tail /var/log/planetgen/planetgen-log`.
- **Behind a reverse proxy** (Cloudflare, a load balancer, nginx in front
  of Apache), set `proxy_fix` in `config.json` first. Without it every
  visitor has the proxy's address, the log shows that address, and the
  jail would ban the proxy: the whole site for everyone.
- Add your own fixed address to `ignoreip` in the jail, and to
  `login_allowlist` in `config.json` if you'd like planetGen itself to
  never lock it either.

## Install (Debian, Ubuntu)

```sh
sudo apt install fail2ban
sudo cp examples/fail2ban/planetgen.filter.conf /etc/fail2ban/filter.d/planetgen.conf
sudo cp examples/fail2ban/planetgen.jail.conf /etc/fail2ban/jail.d/planetgen.conf
# Check the filter against the log (shows how many lines match):
sudo fail2ban-regex /var/log/planetgen/planetgen-log /etc/fail2ban/filter.d/planetgen.conf
sudo fail2ban-client reload
sudo fail2ban-client status planetgen
```

If you set `log_dir` in `config.json`, change `logpath` in the jail to
that folder's `planetgen-log`.

## What the jail does

10 failed sign-ins from one address within 10 minutes ban it for an hour
(`maxretry`, `findtime`, `bantime`). Each later ban of the same address
lasts twice as long as the last, up to a week (`bantime.increment`,
`bantime.factor`, `bantime.maxtime`). Loopback is never banned. The
numbers sit above planetGen's own lockout on purpose: a typo or two never
reaches the firewall, but a script that keeps guessing through the
lockout does.

The log's times are UTC (`2026-10-01T09:43:33Z`); the jail says so with
`logtimezone = UTC`, so it works whatever the server's own time zone is.

To lift a ban: `sudo fail2ban-client set planetgen unbanip <address>`.
To lift planetGen's own lockout as well: Admin › Stats › Login lockouts,
or `python3 -m planetgen.cli.lockouts --ip <address>`.

## Apache-level options

Apache modules can also throttle or filter requests before they reach
planetGen. Both are optional, and neither is set up by `install.sh`:

- **mod_evasive** blocks an address that sends too many requests in a
  short time. It counts every request, not failed sign-ins, so the Galaxy
  Map's many tile requests can trip it; raise its limits well above a
  normal page load if you try it.
- **mod_security** (with the OWASP Core Rule Set) is a web application
  firewall. It guards against attacks planetGen already handles (its
  forms are CSRF-protected and its SQL is parameterized), and its rules
  need tuning to avoid blocking normal use.

planetGen's own lockout plus this fail2ban jail cover password guessing
without either.

## Other platforms

fail2ban is Linux only. On macOS, planetGen's own lockout is
the protection; the activity log has the same lines for any other tool
that reads logs.
