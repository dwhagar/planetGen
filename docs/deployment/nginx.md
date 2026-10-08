# nginx + gunicorn on Linux

nginx serves `/static/` and handles HTTPS. gunicorn runs the Flask app
(every page and the API) as a systemd service. Everything else is the
same as the [Apache setup](apache.md): the checkout at
`/var/lib/planetGen`, `config.json`, `install.sh` and `update.sh`, the
tile cache, the Generate jobs directory, the debug log and the
maintenance timers.

Written for Debian and Ubuntu. For other distributions, see the end of
this page.

## Files

| File | Goes to | Purpose |
|---|---|---|
| [`examples/nginx/planetgen.conf`](../../examples/nginx/planetgen.conf) | `/etc/nginx/sites-available/planetgen` | The site: HTTP to HTTPS redirect, TLS, `/static/` from disk, everything else to gunicorn |
| [`examples/systemd/planetgen-gunicorn.service`](../../examples/systemd/planetgen-gunicorn.service) | `/etc/systemd/system/` | gunicorn: 1 process, 5 threads, on the Unix socket `/run/planetgen/gunicorn.sock`, as `www-data` |
| [`examples/systemd/planetgen-update.service.d/reload-web.conf`](../../examples/systemd/planetgen-update.service.d/reload-web.conf) | `/etc/systemd/system/planetgen-update.service.d/` | Optional: reload gunicorn after the monthly `update.sh` run |

## Install

1. **Database and checkout.** Create the MySQL database and account
   ([`README.md`](README.md#mysql-accounts)), clone the repo to
   `/var/lib/planetGen`, and copy `config.json.example` to `config.json`.
   Fill in `mysql.*` and `secret_key`
   (`python3 -c "import secrets; print(secrets.token_hex(32))"`), and set
   the proxy setting so the app sees the client's address and scheme
   rather than nginx's ([`README.md`](README.md#behind-a-reverse-proxy-proxy_fix)):

       "proxy_fix": {"x_for": 1, "x_proto": 1, "x_host": 0}

2. **`sudo ./install.sh`** from the checkout. It installs the Python
   libraries, migrates the database, fetches the NLTK corpus, and sets
   up the permissions, the tile cache, the jobs directory and the debug
   log. Without Apache installed, its Apache step prints a warning and
   carries on, and it gives the runtime directories to `www-data`, the
   account gunicorn runs as. Ignore its closing Apache instructions.
   Note the admin password it prints in step 2/8. It also offers the
   optional population pass (y/N, default N; see
   [`apache.md`](apache.md#migrating-or-deleting-the-database)).

3. **gunicorn** from apt, so it uses the same system Python and
   libraries `install.sh` set up:

       sudo apt install gunicorn
       gunicorn --version

   The unit runs `gunicorn --pythonpath /var/lib/planetGen/src/html
   wsgi:application`, the same `wsgi.py` Apache's mod_wsgi loads. With
   gunicorn 25.1 or later (Ubuntu 26.04 ships 25.1.0), add
   `--no-control-socket` to `ExecStart`: the control socket tries to
   create `~/.gunicorn` in `www-data`'s home and leaks a thread and a file
   descriptor on every reload in 25.1 to 26.x. Older gunicorn (Ubuntu
   24.04 ships 20.1.0) rejects the option.

       sudo cp examples/systemd/planetgen-gunicorn.service /etc/systemd/system/
       sudo systemctl daemon-reload
       sudo systemctl enable --now planetgen-gunicorn
       systemctl status planetgen-gunicorn

4. **nginx and the certificate.** Install the site and set
   `server_name`. The example points at certbot's certificate paths, so
   get the certificate before enabling the HTTPS `server` block (or
   comment it out until you have one):

       sudo apt install nginx certbot python3-certbot-nginx
       sudo cp examples/nginx/planetgen.conf /etc/nginx/sites-available/planetgen
       sudo ln -s /etc/nginx/sites-available/planetgen /etc/nginx/sites-enabled/
       sudo rm /etc/nginx/sites-enabled/default      # if nothing else uses it
       sudo certbot --nginx -d planetgen.example.com
       sudo nginx -t && sudo systemctl reload nginx

   `http2 on;` needs nginx 1.25.1 or later. On older nginx (Ubuntu 24.04
   ships 1.24), delete that line and write `listen 443 ssl http2;`.

5. **Apache.** If Apache is also installed and was serving the site,
   turn it off so the two don't fight over ports 80 and 443:
   `sudo a2dissite planetgen && sudo systemctl disable --now apache2`.

6. Log in at `https://<your site>/login` and change the admin username
   and password. Then run the [server checklist](../server-checklist.md).

## What the nginx site does

| Apache vhost | nginx site |
|---|---|
| `Alias /static/ .../src/html/static/` | `location /static/ { alias .../src/html/static/; }` |
| `?v=` rule: a year and `immutable`, else `no-cache` | `map $arg_v $planetgen_static_cache` |
| `X-Content-Type-Options: nosniff` on `static/` | `add_header X-Content-Type-Options "nosniff" always;` on `/static/` |
| `AddOutputFilterByType DEFLATE ...` | `gzip on;` with the same types |
| `WSGIDaemonProcess ... processes=1 threads=5` | gunicorn `--workers 1 --worker-class gthread --threads 5` |
| `request-timeout=60` | `proxy_read_timeout 60s` (see below) |
| certbot `--apache` | certbot `--nginx`, plus the HTTP to HTTPS redirect |
| Nothing needed | `X-Forwarded-For`, `X-Forwarded-Proto` headers, read by `proxy_fix` |

The app adds its own security headers (CSP, `X-Frame-Options`,
`Referrer-Policy`, and HSTS on HTTPS). Don't add them in nginx too.

`client_max_body_size 1m`: the API only takes small JSON bodies. Raise it
if a POST fails with 413.

**Timeouts.** mod_wsgi's `request-timeout=60` restarts the daemon when a
request runs past 60 seconds. gunicorn has nothing that exact: with
threads, `--timeout` only catches a worker that stops responding
entirely. nginx gives up after 60 seconds and returns 504, but the thread
keeps working until the request ends.

**Why the socket matters.** `proxy_fix` trusts `X-Forwarded-For` from
whoever connects to gunicorn. The unit binds a Unix socket that only
`www-data` (nginx's group) can open (`--umask 0007`, a 0750 runtime
directory), so only nginx can connect. Don't bind gunicorn to a public
address with `proxy_fix` on.

## Keeping it running

- **After `sudo ./update.sh`:** `sudo systemctl reload
  planetgen-gunicorn`. `update.sh` itself says "reload Apache". For the
  monthly update timer, install the drop-in so it reloads gunicorn for
  you:

      sudo mkdir -p /etc/systemd/system/planetgen-update.service.d
      sudo cp examples/systemd/planetgen-update.service.d/reload-web.conf /etc/systemd/system/planetgen-update.service.d/
      sudo systemctl daemon-reload

- **Reload, not restart.** A reload starts a new worker with the new code
  and leaves a running Generate job alone. `restart` or `stop` ends every
  process in the service, including a Generate job, which the page then
  shows as interrupted. Restarting Apache does the same.
- **Maintenance timers:** unchanged. `examples/maintenance/` is plain
  systemd and doesn't depend on the web server.
- **Logs:** `journalctl -u planetgen-gunicorn` for gunicorn and Python
  errors, `/var/log/nginx/planetgen_*.log` for requests, and the debug
  log as before.
- **More workers:** the rate limits and login lockouts count on Redis by
  default (`ratelimit.storage_uri` empty), so every worker shares them. With
  `memory://` each worker counts on its own ([`api.md`](../api.md#rate-limiting)).

## Checks

The [server checklist](../server-checklist.md) applies, with the gunicorn
commands it lists. And these, specific to the proxy:

    curl -sI "https://HOST/static/style.css?v=1" | grep -i -e cache-control -e nosniff
    # public, max-age=31536000, immutable / nosniff
    curl -sI https://HOST/ | grep -i strict-transport
    # present: proxy_fix is on and the app knows the request was HTTPS
    curl -sI http://HOST/ | head -1
    # 301

## Other distributions

- RHEL, Fedora, Rocky: nginx runs as `nginx` and there is no
  `www-data`. `install.sh` and `update.sh` use apt, so install the
  Python libraries yourself (`pip install ".[api]"` in a virtual
  environment, then point `ExecStart` at its gunicorn), create a
  dedicated account or change `User=`, `Group=` and the ownership steps,
  and add the `nginx` user to that group so it can read `static/` and
  the socket. SELinux also needs the socket labelled for nginx
  (`httpd_var_run_t`) or `setsebool -P httpd_can_network_connect 1` if
  you use a TCP port instead.
- Put the site in `/etc/nginx/conf.d/planetgen.conf` where there is no
  `sites-available`.
