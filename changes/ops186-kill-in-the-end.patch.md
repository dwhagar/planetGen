**Headline:** a worker killed while ending a stage no longer fails the run when a retry is left to finish the end.

### Fixed

- OPS.186: a unit that completed its stage and was killed (`kill -9`, not only the memory killer) in the end of it is retried by RQ, which finishes the end; the run is failed only when no retry is left. Found by Lane 07's galaxy kill test on a second box.
