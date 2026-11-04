**Headline:** update.sh waits up to 90 seconds for the site to answer after a reload, says why each ask failed, and tells a slow start from a server that cannot reach its own address.

### Changed

- **A longer wait (OPS.94).** After the Apache or gunicorn reload, the health check waited 20 seconds. On a real server the first workers were still compiling and importing after new code, each ask timed out at 5 seconds and read as "HTTP 000", and the update opened the site with a warning. `HEALTH_WAIT_SECONDS` now defaults to 90. The wait still ends at the first good answer, so a quick start costs nothing.
- **Every retry line says what curl saw.** "no good answer yet (HTTP 000: the connection was refused: nothing listens on that address and port)", or a timeout, a failed name lookup, a TLS error, instead of a bare code. A line is printed for the first miss, then every 5 seconds or when the reason changes.
- **A server that cannot reach its own address is told apart from a slow start.** When the last answer was no HTTP answer at all and the address is a name, the warning adds what one more ask said with that name pinned to this server's own loopback address: it answered (the server is up, so it cannot reach its own public address: check DNS, the firewall or the router's loopback handling, or set `PLANETGEN_HEALTH_URL`), it did not answer either (the workers had not finished starting, or the server is down), or it answered with an error. This ask only explains the warning; the check still asks the configured address alone (OPS.86), and an address that is an IP number is not asked twice.

generation-output: unchanged
