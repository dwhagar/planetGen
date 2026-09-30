# Windows

planetGen runs natively on Windows with some limits. Three setups are
described here, picked for being maintained, documented, and able to run
unattended as services:

1. **IIS + HttpPlatformHandler**, with IIS starting waitress. For Windows
   Server machines that already run IIS.
2. **Caddy + waitress**, both as Windows services. The simplest to run.
3. **Apache httpd for Windows (Apache Lounge) + waitress** as a service.
   Closest to the Linux [Apache guide](apache.md).

Then [WSL2](#wsl2-the-linux-guides-on-windows), which runs the Linux
guides unchanged.

## Limits on native Windows

- **`install.sh`, `update.sh` and the systemd timers are Linux-only.**
  The steps below do their work by hand.
- **gunicorn does not run on Windows** (it needs `fork`). All three
  setups use waitress, a pure-Python WSGI server that does.
- **The admin Generate page does not work reliably.** Its background
  jobs use POSIX-only process calls: `start_new_session` when starting
  the job runner and each step, `os.killpg` and SIGTERM to cancel, and
  `/proc` or `os.kill(pid, 0)` to check whether a job is still running.
  On Windows, Cancel kills the job runner but leaves the generation step
  running and writing to the database; a finished or crashed job can
  show as still running, which blocks new jobs, or the reverse, letting
  two run at once. An IIS app pool recycle or a service restart also
  ends a running job. Nothing in the web server setup fixes this. The
  code fix is tracked as [`TODO.md`](../TODO.md) item 55.
  **Workaround:** on native Windows, generate from the command line
  (`python generate.py ...` in the checkout, see the
  [README](../../README.md#usage)) and don't use the Generate page. If you
  need the Generate page, use [WSL2](#wsl2-the-linux-guides-on-windows)
  or a Linux VM.
- Everything else works: the pages, the API, the Galaxy Map, admin
  login, and the account and stats pages (the stats page shows no memory
  or load figures, which come from Linux-only sources).

## Which setup

| | 1. IIS + HttpPlatformHandler | 2. Caddy + waitress service | 3. Apache Lounge + waitress service |
|---|---|---|---|
| Best for | Windows Server already running IIS | Anything else; the simplest to run | People who know the Linux Apache guide |
| HTTPS | IIS bindings (certificate store; win-acme for Let's Encrypt) | Automatic | Certificate files (win-acme can write them) |
| Who runs Python | IIS starts and restarts it | A Windows service (WinSW) | A Windows service (WinSW) |
| Extra installs | HttpPlatformHandler, URL Rewrite | Caddy, WinSW | Apache Lounge httpd and its VC++ runtime, WinSW |
| Verified here | No (config written from Microsoft's docs) | Caddyfile checked with `caddy adapt`; waitress tested on Linux | Not run; config follows the Linux vhost |

**Why these three.** Microsoft's own Python-on-IIS guide recommends
HttpPlatformHandler and says wfastcgi is no longer maintained. Caddy
ships native Windows builds with automatic HTTPS and documents running as
a Windows service. Apache Lounge is the Windows build of Apache httpd
that mod_wsgi's documentation points to. **nginx for Windows** was left
out: nginx.org calls it a beta, and it uses only one worker for real
work. **mod_wsgi on Windows** was left out as the app server: it has no
daemon mode on Windows (the app runs inside Apache's own process), its
Windows support is "provisional" in the 6.x releases, installing it
needs Build Tools for Visual Studio, and `jobs.python` must then be set
because `sys.executable` is `httpd.exe`. waitress behind `mod_proxy`
avoids all of that.

## Common steps (all three)

Layout used below (short names, no spaces):

| Path | What |
|---|---|
| `C:\srv\planetGen` | The git checkout (Linux: `/var/lib/planetGen`) |
| `C:\srv\planetgen-venv` | Python virtual environment |
| `C:\ProgramData\planetgen\tiles` | Galaxy Map tile cache (Linux: `/var/cache/planetgen/tiles`) |
| `C:\ProgramData\planetgen\jobs` | Generate jobs (Linux: `/var/lib/planetgen/jobs`) |
| `C:\ProgramData\planetgen\logs` | Debug log, service and web server logs |
| `C:\ProgramData\planetgen\nltk_data` | NLTK `words` corpus (Linux: `/usr/local/share/nltk_data`) |

1. **Database.** Install MySQL 8.4 LTS or MariaDB (current LTS) with its
   Windows installer; both run as a Windows service. Create the database
   and account as in [`README.md`](README.md#mysql-accounts).

2. **Python.** Install Python 3.12 or later from python.org "for all
   users" (so it lands in `C:\Program Files\Python312`, readable by
   service accounts), and Git for Windows. In an elevated PowerShell:

       git clone https://github.com/dwhagar/planetGen.git C:\srv\planetGen
       py -3.12 -m venv C:\srv\planetgen-venv
       C:\srv\planetgen-venv\Scripts\python.exe -m pip install "C:\srv\planetGen[api]" waitress
       C:\srv\planetgen-venv\Scripts\python.exe -m pip uninstall -y planetGen

   The install pulls in the libraries `setup.py` lists. The uninstall
   removes the installed copy of planetGen itself, because the app runs
   from the checkout (as on Linux).

3. **NLTK corpus**, where every account can read it:

       mkdir C:\ProgramData\planetgen\nltk_data
       C:\srv\planetgen-venv\Scripts\python.exe -c "import nltk; nltk.download('words', download_dir=r'C:\ProgramData\planetgen\nltk_data')"

   Each service below sets `NLTK_DATA` to that folder. If you run the
   command-line tools yourself, set it machine-wide too:
   `setx /M NLTK_DATA C:\ProgramData\planetgen\nltk_data`.

4. **config.json.** Copy `config.json.example` to
   `C:\srv\planetGen\config.json` and fill it in. Set the paths
   explicitly: the defaults (`/var/...`) are Linux paths. The changes
   from the example are in
   [`examples/windows/config.json.example`](../../examples/windows/config.json.example):

       "log_file":   "C:/ProgramData/planetgen/logs/planetgen.log",
       "proxy_fix":  { "x_for": 1, "x_proto": 1, "x_host": 0 },
       "tile_cache": { "dir": "C:/ProgramData/planetgen/tiles" },
       "jobs":       { "dir": "C:/ProgramData/planetgen/jobs",
                       "python": "C:/srv/planetgen-venv/Scripts/python.exe" },

   plus `mysql.*` and a `secret_key`
   (`python -c "import secrets; print(secrets.token_hex(32))"`). Forward
   slashes work and avoid JSON escaping. `proxy_fix` makes the app take
   the client's address and scheme from the proxy's `X-Forwarded-For`
   and `X-Forwarded-Proto` headers: rate limits are per client address,
   and the app only sends HSTS when it knows the request was HTTPS (see
   [`README.md`](README.md#behind-a-reverse-proxy-proxy_fix)). For
   option 1 read the note there first.

5. **Schema and first admin login:**

       cd C:\srv\planetGen
       C:\srv\planetgen-venv\Scripts\python.exe src\migrateDb.py

   The first run prints the `admin` password once. Keep it; you change
   it at `/login` later.

6. **Try waitress by hand.** It loads `src\html\wsgi.py`, the file
   Apache's mod_wsgi loads, as the module `wsgi`:

       cd C:\srv\planetGen\src\html
       C:\srv\planetgen-venv\Scripts\waitress-serve.exe --listen=127.0.0.1:8000 --threads=5 --no-clear-untrusted-proxy-headers wsgi:application
       curl.exe http://127.0.0.1:8000/api/health

   `--threads=5` matches Apache's daemon. waitress deletes
   `X-Forwarded-*` headers by default; `--no-clear-untrusted-proxy-headers`
   keeps them for `proxy_fix`. That is safe because waitress listens on
   127.0.0.1 only, so the local proxy is the only thing that can reach
   it. Stop it with Ctrl+C.

7. **Permissions.** The account that runs Python needs read access to
   the checkout, the venv, Python and the corpus, and write access to the
   runtime folders only. As on Linux, it must not be able to change the
   code. With the account from your option as `$acct`:

       $acct = "NT SERVICE\planetgen"          # options 2 and 3
       # $acct = "IIS AppPool\planetgen"       # option 1
       foreach ($d in "tiles","jobs","logs") { mkdir "C:\ProgramData\planetgen\$d" -Force }
       icacls C:\srv\planetGen /grant "${acct}:(OI)(CI)RX"
       icacls C:\srv\planetgen-venv /grant "${acct}:(OI)(CI)RX"
       icacls C:\ProgramData\planetgen\nltk_data /grant "${acct}:(OI)(CI)RX"
       foreach ($d in "tiles","jobs","logs") { icacls "C:\ProgramData\planetgen\$d" /grant "${acct}:(OI)(CI)M" }
       # config.json holds the database password: remove inherited access, keep admins and the app.
       icacls C:\srv\planetGen\config.json /inheritance:r /grant:r "Administrators:F" "SYSTEM:F" "${acct}:R"

   The web server also needs read access to
   `C:\srv\planetGen\src\html\static`. IIS reads it as `IIS_IUSRS`;
   Caddy and Apache as their service account (LocalSystem by default,
   which can already read it).

## Option 1: IIS + HttpPlatformHandler

IIS starts waitress itself on a port it picks (`%HTTP_PLATFORM_PORT%`),
forwards every request to it, and restarts it if it dies.

1. Turn on IIS with these features: Static Content, Static Content
   Compression, Dynamic Content Compression, Request Filtering. Install
   **HttpPlatformHandler v1.2** and **URL Rewrite 2.1** (both from the
   IIS downloads site).
2. Create `C:\srv\planetgen-site` and copy
   [`examples/windows/web.config`](../../examples/windows/web.config) into
   it. This folder is the site's root, not the checkout, so nothing in
   the checkout is reachable except `/static`.
3. In IIS Manager, or in PowerShell:

       Import-Module WebAdministration
       New-WebAppPool planetgen
       Set-ItemProperty IIS:\AppPools\planetgen managedRuntimeVersion ""
       # Keep waitress (and a running request) alive: no idle shutdown, no timed recycle.
       Set-ItemProperty IIS:\AppPools\planetgen processModel.idleTimeout "00:00:00"
       Set-ItemProperty IIS:\AppPools\planetgen recycling.periodicRestart.time "00:00:00"
       New-Website planetgen -PhysicalPath C:\srv\planetgen-site -ApplicationPool planetgen -HostHeader planetgen.example.com -Port 80
       New-WebVirtualDirectory -Site planetgen -Name static -PhysicalPath C:\srv\planetGen\src\html\static

   Then add an HTTPS binding with a certificate (win-acme can get a
   Let's Encrypt one and renew it). The admin pages need HTTPS.

4. Grant `IIS AppPool\planetgen` the permissions from common step 7, and
   `IIS_IUSRS` read on `C:\srv\planetGen\src\html\static`. If IIS says the
   handlers section is locked:
   `%windir%\system32\inetsrv\appcmd unlock config -section:system.webServer/handlers`.

What `web.config` does, compared with the Apache vhost:

- `/static` is the virtual directory, served by IIS's static file handler
  (the Python handler is removed there), with
  `X-Content-Type-Options: nosniff`. Two URL Rewrite outbound rules set
  `Cache-Control` to a year (`immutable`) when the URL has `?v=`, and to
  `no-cache` otherwise.
- Every other URL goes to waitress, started as
  `python -m waitress ... wsgi:application` with `PYTHONPATH` pointing at
  `src\html`. `requestTimeout="00:01:00"` matches `request-timeout=60`.
- 1 MB request body limit, HTTP to HTTPS redirect, compression on.
- `stdoutLogFile` keeps waitress and Python startup errors in
  `C:\ProgramData\planetgen\logs`.

**Check the forwarded headers before trusting them.** This setup has not
been run on a Windows machine. It is not verified whether
HttpPlatformHandler sets or appends `X-Forwarded-For` and
`X-Forwarded-Proto`. If it passes a client's own `X-Forwarded-For`
through unchanged, `proxy_fix` would let a client choose its rate-limit
address. To check: set `"debug": true`, restart the app pool, and from
another machine run
`curl.exe -H "X-Forwarded-For: 192.0.2.1" https://<site>/api/health`.
The last `API request: GET /api/health ... from <address>` line in the
debug log must show that machine's real address. If it shows
`192.0.2.1`, or `127.0.0.1`, set `x_for` to 0 or use the ARR variant
below. Also check that an HTTPS page has a `Strict-Transport-Security`
header. Turn debug off afterwards.

**ARR variant.** Instead of IIS starting Python, run waitress as a
Windows service (option 2, step 1) and let IIS proxy to it with URL
Rewrite and Application Request Routing 3.0:
[`examples/windows/web-arr.config`](../../examples/windows/web-arr.config).
waitress then keeps running across app pool recycles, and a Generate job
is not cut off by one. Enable the proxy at server level, set its timeout
to 60 seconds, turn off "Include TCP port from client IP" (with the port
in `X-Forwarded-For`, every connection would get its own rate-limit
budget), and allow the server variable `HTTP_X_FORWARDED_PROTO` (the
file's comments list the steps).

Updating: see [Updating](#updating-all-options), then recycle the app
pool (`Restart-WebAppPool planetgen`).

## Option 2: Caddy + waitress as a service

1. **waitress as a service.** Download WinSW 2.x (github.com/winsw/winsw,
   `WinSW-x64.exe`), save it as `C:\srv\winsw\planetgen-waitress.exe`, and
   put [`examples/windows/planetgen-waitress.xml`](../../examples/windows/planetgen-waitress.xml)
   next to it. Grant `NT SERVICE\planetgen` the permissions from common
   step 7, then, elevated:

       C:\srv\winsw\planetgen-waitress.exe install
       sc.exe config planetgen obj= "NT SERVICE\planetgen"
       Start-Service planetgen
       curl.exe http://127.0.0.1:8000/api/health

   `NT SERVICE\planetgen` is a virtual account: its own identity, no
   password. The service runs `waitress-serve.exe` from `src\html` with
   the same arguments as common step 6. NSSM works too:

       nssm install planetgen C:\srv\planetgen-venv\Scripts\waitress-serve.exe --listen=127.0.0.1:8000 --threads=5 --no-clear-untrusted-proxy-headers --max-request-body-size=1048576 wsgi:application
       nssm set planetgen AppDirectory C:\srv\planetGen\src\html
       nssm set planetgen AppEnvironmentExtra NLTK_DATA=C:\ProgramData\planetgen\nltk_data PYTHONPATH=C:\srv\planetGen\src\html

2. **Caddy.** Download `caddy.exe` (caddyserver.com/download) to
   `C:\srv\caddy`, with
   [`examples/windows/Caddyfile`](../../examples/windows/Caddyfile) (set
   the site name) and
   [`examples/windows/caddy-service.xml`](../../examples/windows/caddy-service.xml)
   next to it, and WinSW saved as `C:\srv\caddy\caddy-service.exe`:

       C:\srv\caddy\caddy.exe validate --config C:\srv\caddy\Caddyfile
       C:\srv\caddy\caddy-service.exe install
       Start-Service caddy

   Open ports 80 and 443 in Windows Firewall. With the name pointing at
   this machine, Caddy gets the certificate by itself.

The Caddyfile is the [Linux one](caddy.md) with Windows paths and a TCP
upstream: `/static/` from the checkout with nosniff and the `?v=` cache
rule, compression, a 1 MB body limit, 60 seconds to the first response
byte, and everything else to waitress on 127.0.0.1:8000. Caddy sets
`X-Forwarded-For` and `X-Forwarded-Proto` itself.

Updating: see [Updating](#updating-all-options), then
`Restart-Service planetgen`.

## Option 3: Apache (Apache Lounge) + waitress as a service

1. waitress as a service: option 2, step 1.
2. Install the Apache Lounge build of httpd 2.4 (apachelounge.com) to
   `C:\Apache24`, and the Visual C++ redistributable it names. In
   `C:\Apache24\conf\httpd.conf`, uncomment the `LoadModule` lines for
   `headers`, `deflate`, `proxy`, `proxy_http`, `ssl`, `socache_shmcb`
   and `rewrite`, and add at the end:

       Include conf/extra/httpd-planetgen.conf

3. Copy [`examples/windows/httpd-planetgen.conf`](../../examples/windows/httpd-planetgen.conf)
   to `C:\Apache24\conf\extra\`, set `ServerName` and the certificate
   paths, then:

       C:\Apache24\bin\httpd.exe -t
       C:\Apache24\bin\httpd.exe -k install
       Start-Service Apache2.4

The file is the Linux vhost with two changes: `mod_proxy_http` to
waitress in place of `WSGIScriptAlias`, and an HTTPS virtual host. The
`Alias /static/`, the `<Directory>` block with nosniff and the `?v=`
`<If>`/`<Else>` rule, and `DEFLATE` are the same. `ProxyPass ...
timeout=60` stands in for `request-timeout=60`. `mod_proxy_http` appends
`X-Forwarded-For` itself, and the file sets `X-Forwarded-Proto: https`.

## Updating (all options)

`update.sh` doesn't run on Windows. By hand, elevated:

    cd C:\srv\planetGen
    git pull
    C:\srv\planetgen-venv\Scripts\python.exe -m pip install --upgrade "C:\srv\planetGen[api]" waitress
    C:\srv\planetgen-venv\Scripts\python.exe -m pip uninstall -y planetGen
    C:\srv\planetgen-venv\Scripts\python.exe src\migrateDb.py

then restart the app (per option above).
`python src\migrateDb.py --status` shows whether a migration is pending.
New files inherit the checkout's permissions, so the service account
still only has read access.

## Monthly orbit update (Task Scheduler)

The counterpart of `planetgen-orbits@.timer`. Copy
[`examples/windows/planetgen-orbits.cmd`](../../examples/windows/planetgen-orbits.cmd)
to `C:\srv\`, then create one task per database:

    schtasks /Create /TN "planetGen orbits (planetgen)" /SC MONTHLY /D 1 /ST 03:30 `
        /TR "C:\srv\planetgen-orbits.cmd planetgen" /RU "NT AUTHORITY\SYSTEM"

It reads the database credentials from `config.json`. In Task Scheduler,
tick "Run task as soon as possible after a scheduled start is missed"
(the equivalent of `Persistent=true`). The Linux update timer has no
Windows counterpart: update by hand as above.

## Debug log

With `"debug": true`, the log goes to `log_file`
(`C:/ProgramData/planetgen/logs/planetgen.log` above). Nothing rotates it
on Windows, and it grows fast (megabytes per sector). Turn debug off when
you are done, or delete the file while the service is stopped.

## WSL2 (the Linux guides on Windows)

WSL2 runs a real Linux kernel, so the Linux guides ([Apache](apache.md),
[nginx](nginx.md) or [Caddy](caddy.md)) work unchanged, including the
Generate page and the systemd timers.

1. `wsl --install -d Ubuntu`, then enable systemd in Ubuntu's
   `/etc/wsl.conf` (current Ubuntu images already have it):

       [boot]
       systemd=true

   and run `wsl --shutdown` from Windows to apply it.
2. Follow a Linux guide inside Ubuntu. MySQL can run inside WSL or on
   Windows.
3. To reach it from other machines: on Windows 11 22H2 or later, set
   `networkingMode=mirrored` under `[wsl2]` in `%UserProfile%\.wslconfig`
   so WSL shares the Windows host's addresses; otherwise forward ports
   with `netsh interface portproxy`. Also allow inbound traffic through
   the Hyper-V firewall (Microsoft's WSL networking page has the
   `New-NetFirewallHyperVRule` command).

WSL2 starts when a user starts it. There is no supported way yet to start
it at boot with nobody logged in, so for an unattended server use a Linux
VM (Hyper-V) instead.
