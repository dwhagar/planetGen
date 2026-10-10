# Running CI on your own computers (optional)

The CI workflow runs only when started by hand (Actions > CI > Run
workflow). When it runs, it uses GitHub's machines by default (`ubuntu-latest`,
`macos-latest`), and nothing in this guide is needed for
it to work. Self-hosted runners are **off by default**. This guide is for
turning them back on later: it covers what each machine needs, how to
register it, and the security settings to apply first.

## What runs where

Self-hosted runners are used only when the repository variable
`SELF_HOSTED` is set to `true` (Settings > Secrets and variables > Actions
> Variables). With it unset, or set to anything else, every job runs on
GitHub's machines and the `RUNNER_*` variables below are ignored, so a
leftover variable can't leave jobs queued for a computer that no longer
exists.

With `SELF_HOSTED` set to `true`, each job picks its runner from the
variable in the table. The value is JSON: either a list of labels, such as
`["self-hosted", "Linux"]`, or a GitHub-hosted runner name in quotes, such
as `"ubuntu-latest"`. If the variable is unset, the job uses the
self-hosted default in the table.

| Variable | Jobs | Default when `SELF_HOSTED` is `true` | Needs |
|---|---|---|---|
| `RUNNER_LINUX` | tests (3 database engines), browser checks, dependency audit, deep fuzz, release note, version stamp | `["self-hosted", "Linux"]` | Linux with Docker |
| `RUNNER_LINUX_INSTALLERS` | `install.sh` and `update.sh` on Linux against a live database | `"ubuntu-latest"` | Linux with Docker; see [The installer jobs](#the-installer-jobs) |
| `RUNNER_MACOS_INSTALLERS` | `install.sh` on macOS | `"macos-latest"` | see [macOS](#macos) |

A job runs on any online runner that has **all** of its labels. You can
have several Linux runners: jobs spread across them, and a runner takes
one job at a time. To go back to GitHub's machines, delete `SELF_HOSTED`
or set it to `false`.

A pull request from a fork always runs on GitHub's machines, whatever the
variables say (see [Security](#security)).

## Security

This repository is **public**. Anyone can fork it and open a pull request,
and a workflow runs whatever code that pull request contains. On a
self-hosted runner, that code would run on your computer. Set these up
before you register a runner:

1. **Fork pull requests.** The workflows already send pull requests from
   forks to GitHub-hosted runners. As a second guard, go to Settings >
   Actions > General > "Approval for running fork pull request workflows
   from contributors" and choose **Require approval for all external
   contributors**. Nothing from a fork then runs until you approve it.
2. **Use a dedicated account.** Run the runner as its own user (for
   example `gha`), never your everyday login or an administrator. On
   Linux, membership of the `docker` group gives that user root-equivalent
   access to the machine, so a spare machine or a VM is better than your
   desktop.
3. **Keep secrets off the machine.** Don't leave SSH keys, cloud
   credentials, browser profiles or a `config.json` with real database
   passwords in that user's home directory. The jobs need no secrets.
4. **Leave workflow permissions as they are** (Settings > Actions >
   General > Workflow permissions). Only `stamp-version.yml` writes to the
   repository (it commits the release), and it asks for that permission
   itself.
5. **Check which repositories a runner serves.** A runner added under this
   repository's Settings > Actions > Runners only takes jobs from this
   repository.

## Linux runner

Tested on Ubuntu 22.04 and 24.04 and Debian 12. Any 64-bit distribution
with Docker works, but see "setup-python can't find 3.9" under
[Troubleshooting](#troubleshooting) for distributions other than
Ubuntu.

1. **Create the runner user and install the basics** (as an administrator):

   ```sh
   sudo adduser --disabled-password --gecos "" gha
   sudo apt-get update
   sudo apt-get install -y git curl tar unzip docker.io python3-venv
   sudo usermod -aG docker gha
   ```

2. **Register the runner.** In GitHub, go to Settings > Actions > Runners
   > New self-hosted runner > Linux, and run the commands it shows as
   `gha` (`sudo -iu gha`). Use a folder of its own, such as
   `~/actions-runner`. When `./config.sh` asks:
   - name: something you'll recognise, such as `linux-box-1`;
   - labels: press Enter. The default labels are `self-hosted`, `Linux`
     and `X64` (or `ARM64`), which is what the workflows ask for;
   - work folder: press Enter (`_work`).

   The token on that page is only valid for an hour. If it expires, reload
   the page for a new one.

3. **Install it as a service**, so it starts at boot (as an administrator,
   from the runner folder):

   ```sh
   cd /home/gha/actions-runner
   sudo ./svc.sh install gha
   sudo ./svc.sh start
   ```

   It should now show as **Idle** under Settings > Actions > Runners.

4. **Chromium's system libraries** for the browser checks (once, as an
   administrator). The job installs the browser itself, but it needs root
   for the system libraries and skips that step when `sudo` asks for a
   password:

   ```sh
   python3 -m venv /tmp/pw && /tmp/pw/bin/pip install playwright
   sudo /tmp/pw/bin/python -m playwright install-deps chromium
   rm -rf /tmp/pw
   ```

The jobs set up everything else themselves:

- **Python**: `actions/setup-python` downloads 3.9 and 3.12 into the
  runner's tool cache (`_work/_tool`) the first time they are needed.
- **The database**: each test job starts its own MySQL or MariaDB
  container on a free port, so jobs don't collide. A MySQL server you
  already run on port 3306 is fine too.
- **Packages**: each job installs into its own virtualenv, which is
  deleted when the job ends.
- **NLTK's `words` corpus**: downloaded once to `~/nltk_data`.

## macOS

No job needs a Mac by default. The only macOS job, the installer check,
stays on GitHub's machines (see the next section). A Mac runner is set up
the same way as Linux, using the macOS download on the New self-hosted
runner page and `./svc.sh install` without a user name. It has no Docker
service containers, so it can't take the Linux test jobs.

## The installer jobs

`linux-update` and `macos-installers` run
`install.sh` and `update.sh` with sudo (`linux-update` against a MySQL service container, through a
fresh install, an update with nothing new, a migration, a database newer
than the code, an unreachable server and a failed migration). They copy
planetGen into system folders, install packages, register launchd
daemons, and start a web server on port 8000. On GitHub's machines all of that is thrown away
after the job. On your own computer it stays, and the next run trips over
it. Keep their variables unset unless you have a machine you're happy to
wipe, such as a VM you restore to a snapshot after each run.

## Several runners, and adding a machine

- **A new machine** (your third one, say): follow the section for its OS.
  It joins the pool as soon as it shows Idle; nothing in the workflows
  changes.
- **Several runners on one machine**: register each in its own folder
  (`~/actions-runner-1`, `~/actions-runner-2`) with its own name, and
  install each as its own service. Jobs on one machine don't collide: each
  database container gets its own port, and each job has its own
  virtualenv and work folder. The test job uses every core (`pytest -n
  auto`), so two test jobs on one machine each run at about half speed.
- **Different labels**: if a runner has extra or different labels, set
  `RUNNER_LINUX` to the list of labels to require,
  for example `["self-hosted", "Linux", "X64", "fast"]`.

## Upkeep

- The runner updates itself.
- `actions/checkout` cleans the work folder at the start of every job.
  Docker images and stopped containers do pile up, so run `docker system
  prune -af` now and then (or weekly from cron).
- Runners stay registered (persistent). GitHub also supports `--ephemeral`
  runners that take one job and then unregister, which isolates jobs
  further, but something then has to re-register them after every job.
  That isn't needed for this repository.

## Troubleshooting

| What you see | Cause and fix |
|---|---|
| A job sits in "Queued" forever | No online runner has all of its labels. Compare the job's labels (shown on the job page) with Settings > Actions > Runners. Fix the runner's labels, set the matching `RUNNER_*` variable, or set it to `"ubuntu-latest"` to use GitHub's machines. |
| Runner shows Offline | The service isn't running. Run `sudo ./svc.sh status`, then `start`. |
| `docker: command not found`, or `permission denied` on `/var/run/docker.sock` | Install Docker and add the runner user to the `docker` group, then restart the service (`sudo ./svc.sh stop && sudo ./svc.sh start`). |
| A service container never becomes healthy | Docker can't pull the image (no internet or a proxy), or the machine is out of memory or disk. Run `docker pull mysql:8.0` by hand as the runner user to see the error. |
| setup-python can't find 3.9 or 3.12 | Its downloads are built for Ubuntu. On other distributions, install those versions into the runner's tool cache yourself, or run this runner in an Ubuntu VM. |
| Hundreds of jobs queued, new ones wait for hours | Runs queued before the runners came online are worked through oldest first. Cancel the stale ones: Actions tab > filter "is:queued" > open each run > Cancel workflow run. |
| Old pull requests' "Release note" check fails saying the PR "changes the version" | A run that waited in the queue was comparing against today's main. Fixed: the check now compares against the base the run started from. Re-run it, or ignore it on merged PRs. |
| The browser job fails with missing `.so` libraries | Chromium's system libraries aren't installed. See step 4 under [Linux runner](#linux-runner). |
| Disk full | Run `docker system prune -af` and empty `_work/_tool` of Python versions you no longer need. |
| Port 8000 already in use on an installer job | That job ran on your own machine before and its web server is still running. Stop it and keep the installer jobs on GitHub's machines. |
