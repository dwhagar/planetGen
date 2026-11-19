**Headline:** A job whose arguments are bigger than the control socket takes is handed them by a file, and a job that cannot be started holds its run (OPS.283).

### Fixed

- OPS.283 (bug): the arguments of a unit over 64 KB go to the job in a file in the state folder (`job-arguments`) instead of down the control socket, which refused a message over about 213 KB. The template reads the file and deletes it, the supervisor deletes it when the hand-over fails, and a file left by a supervisor that ended is swept after a day. The template reads only files of that folder.
- OPS.283 (bug): a job that cannot be started (its process was never forked) left its row `failed` while its run went on `running` with a unit that would never end. The run now holds with the error text, the unit among the stage's held ones, as for a unit that failed after it ran. A page step fails the steps after it and an opening holds its run.

generation-output: unchanged
generation-affects: none
