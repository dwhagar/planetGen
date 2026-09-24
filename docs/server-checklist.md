# Server checklist: confirm the site loads

Run this after deploying to confirm the fixes for the "site won't load"
outages are live: CGI pages logging `http/client.py` timeouts,
"Truncated or oversized response headers", and Apache being OOM-killed
after Galaxy Map traffic.

Replace `HOST` with your site's hostname, `DB_NAME` with your database,
and `/var/lib/planetGen` with your checkout if it lives elsewhere.

## 1. Is the new code deployed?

    cd /var/lib/planetGen && git log -1 --oneline
    grep __version__ src/stellarObjects/_version.py

Pass: version `5.47.0` or later. That release carries the tile-based Galaxy
Map (the OOM fix) and schema v26.

## 2. Has the database actually been migrated?

Restarting Apache does not apply schema migrations. Only `migrateDb.py`
(or `update.sh`/`install.sh`, which call it) does.

    sudo systemctl reload apache2
    curl -s https://HOST/api/health

Pass: `"schema_version": 26, "schema_current": true`.
If not: run `python3 src/migrateDb.py` with the same `PLANETGEN_MYSQL_*`
settings the site uses, then check again.

If you use more than one database, check each one. `migrateDb.py` only
migrates the one it is pointed at:

    curl -s "https://HOST/api/health?db=OTHER_DB_NAME"
    python3 src/migrateDb.py --mysql-database OTHER_DB_NAME

## 3. Are the spatial indexes there?

    mysql DB_NAME -e "SELECT table_name, index_name FROM information_schema.statistics
      WHERE table_schema = DATABASE() AND index_name LIKE 'idx_%_center' AND seq_in_index = 1"

Pass: five rows: `sectors` (v25) plus `nebulae`, `asteroid_fields`,
`black_holes` and `neutron_stars` (v26).

## 4. Is the Apache config current?

    grep -n "WSGIDaemonProcess planetgen-api" /etc/apache2/sites-available/planetgen.conf

Pass: the line includes `request-timeout=60` (see
`examples/apache/planetgen.conf.example`). Without it, one runaway request
can hold an API thread forever.

## 5. Is the tile cache writable?

The Galaxy Map caches tiles on disk, in `/var/cache/planetgen/tiles` by
default (see `tile_cache` in [`config.md`](config.md)).

    ls -ld /var/cache/planetgen/tiles

Pass: the folder exists and is owned by (or writable by) Apache's user,
usually `www-data`. It fills up after you open the Galaxy Map.

## 6. Time the Galaxy Map's API calls

    time curl -s -o /dev/null "https://HOST/api/galaxy/stamp"
    time curl -s -o /dev/null "https://HOST/api/galaxy/tiles?tiles=0/0/0/0"

Pass: each finishes in a few seconds at most. CGI pages give up at 30
seconds, which is when they log the `http/client.py` traceback.

## 7. Reproduce real use while watching the logs

    sudo tail -f /var/log/apache2/planetgen_error.log | grep -E "Timeout|Truncated|http/client|MemoryError"

In a second terminal, watch the API's memory:

    watch -n 2 "ps -o rss,cmd -C apache2 | sort -n | tail -3"

Open the Galaxy Map, click empty space well away from the core so the view
re-centers there, then zoom and rotate for a minute. While the map updates,
open another page (for example Browse) in a second tab.

Pass: nothing new in the error log, memory stays flat instead of climbing
by hundreds of MB, and the other page loads normally.
