# Windows

planetGen is not supported or tested on Windows. It targets Linux and
macOS, and nothing in the repository (installers, update scripts, service
files, scheduled-task scripts) is written for Windows. This page only
describes the simplest way a typical Windows machine can still run it.

## The simplest setup: WSL2 with Ubuntu

Install [WSL2](https://learn.microsoft.com/windows/wsl/install) with an
Ubuntu distribution, then work entirely inside it as you would on a Linux
machine:

1. Open the Ubuntu shell and follow [INSTALL.md](../INSTALL.md) from the
   top: Python 3.9 or later, MySQL 8.0.16+ or MariaDB 10.4+, Redis, and a
   git checkout of the repository.
2. For the website, set up Apache with mod_wsgi as in
   [deployment/apache.md](deployment/apache.md), or another Linux option
   from [deployment/README.md](deployment/README.md).
3. Keep the checkout on the Linux filesystem (somewhere under your Ubuntu
   home folder or `/var/lib/planetGen`), not under `/mnt/c`, which is slower.

WSL2 may stop an idle distribution, and with it MySQL/MariaDB and Redis.
The command line and the website need both running whenever you use them.

## Running from the checkout

The code runs only from the git checkout, installed as an editable package:

```bash
pip install -e .
```

For a quick look at the website without a web server, start the
development server from the checkout:

```bash
python src/html/wsgi.py
```

Before that, `config.json` (copied from `config.json.example`) must point
at a running MySQL or MariaDB and a running Redis; see
[config.md](config.md). The development server is for trying things out,
not for serving a site.

## Native Windows

There is no native Windows install. Anyone who wants one has to do the
work themselves: obtaining Python, MySQL or MariaDB and Redis for Windows,
running a WSGI server, and arranging services and scheduled maintenance.
No installer, service files or scheduled-task scripts are provided.
