**Headline:** When a long sky job loses its MySQL connection, the log now says why: a server restart, an idle timeout, an abort by the server, or something closing the socket.

### Changed

- A sky unit (`run_galaxy._sky_unit`) reads the server's uptime, timeouts and aborted-client count when it opens its connection (`db/connection_loss.py`). When the connection is lost it reads them again on the new one and logs one sentence, for example `MySQL connection lost: MySQL was restarted: its uptime is 40 s now and was 5,000 s 1,800 s ago`, with the logs to look in. The other verdicts are an idle timeout (no batch finished for as long as `wait_timeout`), the server aborting the client (`Aborted_clients` went up: a read or write timeout or a packet over `max_allowed_packet`), and the same server with something between the job and MySQL closing the socket.
- Nothing about the retry changes: a unit still goes on from its last committed batch on a new connection, up to 3 times. An account that cannot read the server status gets a line that says it cannot tell. DB.115 stays open until a real loss logs a verdict.

generation-output: unchanged
