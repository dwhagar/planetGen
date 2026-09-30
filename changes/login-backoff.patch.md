### Security
- **Failed logins are counted per username, not only per address.** After
  10 failures for one username, from any mix of addresses, each further
  failure locks that username for twice as long as the last (1 s, 2 s,
  4 s, ... up to 15 minutes). A locked login is refused with a 429 and
  `Retry-After` before the password is checked, and the login page says
  how long to wait. Unknown usernames are counted the same way, so the
  lock doesn't reveal which usernames exist; a successful login clears
  the count (`src/html/api/loginbackoff.py`).
