**Headline:** INSTALL.md walks through giving Redis a user and a password, for a new install and for one already running without.

### Changed

- `INSTALL.md` (DOC.109): a new section, "Give Redis a user and a password", with the steps in order: choose a password, copy `examples/redis/planetgen-acl.conf` and put the password in its `user planetgen` line, include it at the end of `redis.conf` (which also switches the `default` user off), put the password in `redis.url`, restart Redis and the web server, check with `redis-cli` and `update.sh`, and how to go back. It says what the `update.sh` warning means and that it goes away once the connection logs in as the `planetgen` user.
- `docs/deployment/README.md`: the Redis settings entry no longer says a password can wait; it points to those steps.

- `INSTALL.md`: "cancelled" is spelled "canceled" (the spelling test was red on main for that word).

generation-output: unchanged
