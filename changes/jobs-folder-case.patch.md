### Fixed

- The Generate page's jobs folder defaults to `/var/lib/planetGen/jobs`, in the checkout's own folder, instead of `/var/lib/planetgen/jobs`, which differed from it only by case and left two folders; `update.sh` moves jobs from the old folder into the new one (a running job stays put) and removes the old folder once empty. A `jobs.dir` set in `config.json` is left alone (OPS.19).
