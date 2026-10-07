# Caddy + gunicorn on Linux

Caddy in place of nginx or Apache. It gets and renews the TLS certificate
itself, redirects HTTP to HTTPS without extra config, and the whole site
fits in one short file. gunicorn runs the app exactly as in the
[nginx guide](nginx.md); only the front end changes.

## Files

| File | Goes to |
|---|---|
| [`examples/caddy/Caddyfile`](../../examples/caddy/Caddyfile) | `/etc/caddy/Caddyfile` |
| [`examples/systemd/planetgen-gunicorn.service`](../../examples/systemd/planetgen-gunicorn.service) | `/etc/systemd/system/` |

## Install

1. Follow the [nginx guide](nginx.md#install), steps 1 to 3: the
   database, the checkout, `config.json` with
   `"proxy_fix": {"x_for": 1, "x_proto": 1, "x_host": 0}`,
   `sudo ./install.sh`, and the `planetgen-gunicorn` service. Skip nginx.

2. Install Caddy from its own apt repository (the steps are at
   caddyserver.com/docs/install). It comes with a `caddy` systemd service
   that runs as the `caddy` user.

3. Let Caddy read `static/` and connect to gunicorn's socket. Both belong
   to the `www-data` group (the code tree is `root:www-data`, mode
   750/640):

       sudo usermod -aG www-data caddy

4. Install the Caddyfile, set the site name, check it, and restart
   Caddy (a restart, not a reload, so the new group membership applies):

       sudo cp examples/caddy/Caddyfile /etc/caddy/Caddyfile
       sudo mkdir -p /var/log/caddy && sudo chown caddy:caddy /var/log/caddy
       sudo -u caddy caddy validate --config /etc/caddy/Caddyfile
       sudo systemctl restart caddy

   The name must resolve to this server, and ports 80 and 443 must be
   reachable, or Caddy can't get a certificate. `journalctl -u caddy`
   shows how that went.

5. Log in at `https://<your site>/login` and change the admin username
   and password. Then run the [server checklist](../server-checklist.md).

## What the Caddyfile does

| Apache vhost | Caddyfile |
|---|---|
| `Alias /static/ .../src/html/static/` | `handle_path /static/* { root * .../src/html/static; file_server }` |
| `?v=` Cache-Control rule | `@versioned query v=*` and `@unversioned not query v=*`, each with a `header` |
| `X-Content-Type-Options: nosniff` on `static/` | `header X-Content-Type-Options "nosniff"` in the static block |
| `mod_deflate` | `encode zstd gzip` |
| `request-timeout=60` | `response_header_timeout 60s` |
| certbot | Automatic |

Caddy sets `X-Forwarded-For`, `X-Forwarded-Proto` and `X-Forwarded-Host`
on every proxied request and drops any a client sent, so `proxy_fix`
gets the real client address and scheme. The app adds its own security
headers, including HSTS on HTTPS.

The timeout caveat from the nginx guide applies: Caddy gives up after 60
seconds, but the gunicorn thread keeps working until the request ends.

## Keeping it running

As in the [nginx guide](nginx.md#keeping-it-running):
`sudo systemctl reload planetgen-gunicorn` after `update.sh` (or install
the drop-in for the monthly timer), reload rather than restart while a
Generate job runs, and the maintenance timers are unchanged. Caddy's own
config reload is `sudo systemctl reload caddy`.

## Checks

    curl -sI "https://HOST/static/style.css?v=1" | grep -i -e cache-control -e nosniff
    curl -sI https://HOST/ | grep -i strict-transport
    curl -s https://HOST/api/health

Then the [server checklist](../server-checklist.md).
