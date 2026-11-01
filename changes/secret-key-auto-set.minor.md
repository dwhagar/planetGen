**Headline:** A missing secret key is now generated and kept, instead of a new random one for every process.

### Changed
- With no `secret_key` in `config.json` or `PLANETGEN_SECRET_KEY`, the update (and, as a fallback, the app at start) says a key is missing, makes one, and keeps it in the state folder (`secret_key`, mode 0600). Every worker and restart reads the same key, so forms left open across a restart keep working. The key is never printed or logged. An explicit setting still wins.
