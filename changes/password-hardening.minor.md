### Added

- Admin sign-in hardening (SEC.22 to SEC.25): a wrong current password on
  the account page now counts as a failed login and is limited like one;
  a successful login gives the browser a 90-day trusted-device cookie
  that keeps it past its username's lockout (the per-address lockout
  still applies); new passwords are checked against a bundled list of
  about 47,000 common and breached passwords and may not be just the
  username or the site's name with a few characters added; passwords are
  hashed with PBKDF2-SHA256 at 600,000 rounds and older hashes are
  upgraded at the next login. `src/loginLockouts.py --forget-devices
  <user>` revokes an admin's trusted devices. Control schema v3 (rerun
  `update.sh`).
- A fail2ban filter and jail for failed admin sign-ins
  (`examples/fail2ban/`, guide in `docs/deployment/fail2ban.md`; SEC.27).
- Optional two-factor sign-in for admins (SEC.26): turn it on from the
  account page by scanning a QR code with any authenticator app; sign-in
  then asks for the 6-digit code after the password (or one of ten
  single-use recovery codes). Wrong codes count as failed logins and are
  logged as `AUTH totp.failed` (matched by the fail2ban filter). API keys
  are unaffected. `src/loginLockouts.py --reset-two-factor <user>` turns
  it off for a lost phone. Control schema v4 (rerun `update.sh`).
