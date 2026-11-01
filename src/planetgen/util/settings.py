# planetgen/util/settings.py

"""
The settings model (ADM.42)
===========================

One Pydantic model, `Settings`, describes every `config.json` option: its
type, default, help text, unit, whether it is a secret, whether changing it
needs a restart, which environment variable overrides it and whether the
web may edit it. The loader (`get_settings`), the generated tables in
`docs/reference/config.md`, `config.json.example` and `config.schema.json`
(`planetgen.cli.config`) all read the model, and tests fail when they
drift or when code reads a `PLANETGEN_*` variable the model doesn't name.

One file, one model (Boss, 2026-10-16: "1 file for all program configuration
settings"): `config.json` in the install folder, which the installers give to the
web server's user, holds every option, and the configuration page and the command
line (`planetgen.util.settings_writer`) edit it. Load order, lowest to highest: the
model's defaults, `config.json`, the `PLANETGEN_*` environment variables, then an
explicit argument at the call site.

Field metadata sits in `json_schema_extra` under `x-` keys, written with
`opt` below:

- `x-env`: the environment variable that overrides it
- `x-env-empty`: an environment variable that is set but empty still
  overrides (a blank MySQL password), instead of counting as unset
- `x-restart`: a change only applies after a restart
- `x-secret`: never shown in a page or log
- `x-editable`: `false` means the web may not change it (the database, secrets
  that open the host, `jobs.python`); only the command line, run on the server, can
- `x-min_role`: `admin` (the default) or `owner`: who may save it from the web
- `x-multiline`: a text option the form edits in a text area (robots.txt text, security.txt)
- `x-unit`: the form's suffix
- `x-label`, `x-summary`: the option's name in plain words and one sentence on what it does (the long
  `description` is the help under "Details"); both are required, and a test holds them to the page's wording rules
- `x-group`: the purpose-based group the settings page lists it under (`planetgen.util.settings_view.GROUPS`)
- `x-control`: the control the page draws when the type alone doesn't say (`slider`, `duration`, `rate`, `color`,
  `asset`, `language`); empty means the control the type gives (switch, choice, number, text, chips, text area)
- `x-advanced`: shown only under "Advanced" in its group (search always reaches it)
- `x-choice-labels`: `{value: words}` for a choice option whose values are not words an admin reads
- `x-slider`: `(low, high)`, the range a `slider` control's logarithmic stops span (0 is added when the option allows it)

The log options (`debug`, `log_file`, `log_dir`, `log_rotation`) are read
by `planetgen.util.logpaths`, which needs only the standard library (the
installers run it as root before any virtual environment exists); their
defaults come from there, and a test checks the two agree.

Annotations here are `Optional[...]`/`List[...]`, as in most of the package.
"""

import copy
import json
import os
import re
import threading
from typing import Any, List, Literal, get_args, get_origin
from urllib.parse import urlsplit

from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator

from planetgen.util import logpaths

CONFIG_VERSION = 1
"""int: The `config_version` this code writes and reads. A file without one
is version 1; a newer one is refused."""

ENV_PREFIX = "PLANETGEN_"

_FALSE_STRINGS = ("", "0", "false", "no", "off")


def opt(default, description, env=None, restart=False, secret=False, editable=True, role="admin",
        unit=None, env_empty=False, multiline=False, *, label, summary, group, control=None, advanced=False,
        choice_labels=None, slider=None, **constraints):
    """A model field with its `x-` metadata (see the module docstring)."""
    extra = {"x-restart": restart, "x-secret": secret, "x-editable": editable, "x-min_role": role,
             "x-label": label, "x-summary": summary, "x-group": group, "x-control": control,
             "x-advanced": advanced, "x-multiline": multiline, "x-choice-labels": choice_labels, "x-slider": slider}
    if env:
        extra["x-env"] = env
    if env_empty:
        extra["x-env-empty"] = True
    if unit:
        extra["x-unit"] = unit
    return Field(default, description=description, json_schema_extra=extra, **constraints)


def _loose_bool(value):
    """`false`, `0`, `no`, `off` and empty text mean off, any other text on
    (what `PLANETGEN_DEBUG` always did)."""
    if isinstance(value, str):
        return value.strip().lower() not in _FALSE_STRINGS
    return value


class _Section(BaseModel):
    model_config = ConfigDict(extra="ignore", validate_default=True)


class MySQL(_Section):
    host: str = opt("127.0.0.1", "The MySQL server's host name or address.", env="PLANETGEN_MYSQL_HOST",
                    restart=True, editable=False, group="database", env_empty=True,
                    label="Database host",
                    summary="The address of the database server.")
    port: int = opt(3306, "The server's port.", env="PLANETGEN_MYSQL_PORT", restart=True, editable=False,
                    group="database", ge=1, le=65535, env_empty=True,
                    unit="port number", label="Database port",
                    summary="The port the database server listens on.")
    user: str = opt("planetgen", "The account every entry point connects with; give it whatever grants the most "
                    "demanding caller needs.", env="PLANETGEN_MYSQL_USER", restart=True, editable=False,
                    group="database", env_empty=True,
                    label="Database user",
                    summary="The account every part of planetGen connects with.")
    password: str = opt("", "That account's password.", env="PLANETGEN_MYSQL_PASSWORD", restart=True, secret=True,
                        editable=False, group="database", env_empty=True,
                        label="Database password",
                        summary="That account's password; it is never shown.")
    database: str = opt("planetgen", "The galaxy database the pages show and the generators fill.",
                        env="PLANETGEN_MYSQL_DATABASE", restart=True, editable=False, group="database",
                        env_empty=True,
                        label="Galaxy database",
                        summary="The database the pages show and the generators fill.")
    database_prefix: str = opt("planetgen", "The schema-name prefix the database picker lists, for a deployment "
                               "with several galaxies on one server.", env="PLANETGEN_MYSQL_DATABASE_PREFIX",
                               restart=True, editable=False, group="database",
                               label="Galaxy database prefix",
                               summary="The name prefix the database picker lists, for several "
                                       "galaxies on one server.")
    statement_timeout_seconds: float = opt(
        10, "The longest any statement on the web interface's and API's read-only connections may run (0 turns "
        "the limit off). A query it stops gives a 504 page, or a QUERY_TIMEOUT error from the API.",
        unit="seconds", group="housekeeping", ge=0,
        label="Longest page query",
        summary="How long a page or API query may run before it is stopped; 0 means no limit.",
        control="duration")


class Redis(_Section):
    url: str = opt("redis://127.0.0.1:6379/0", "The Redis server the work queue, the rate limits and the login "
                   "lockouts use.", env="PLANETGEN_REDIS_URL", restart=True, editable=False, group="queue",
                   label="Redis address",
                   summary="The Redis server used by the work queue, the rate limits and the login "
                           "lockouts.")


class RatePages(_Section):
    search: str = opt("30 per minute", "Limit per client address on /search and /galaxy/locate. Empty turns it off.",
                      restart=True, group="protection",
                      label="Search limit",
                      summary="How often one visitor may search; empty turns the limit off.",
                      control="rate")
    galaxy: str = opt("60 per minute", "Limit on /galaxy, the Galaxy Map page. Empty turns it off.", restart=True,
                      group="protection",
                      label="Galaxy Map page limit",
                      summary="How often one visitor may open the Galaxy Map page; empty turns the "
                              "limit off.",
                      control="rate")
    galaxy_tiles: str = opt("600 per minute", "Limit on /galaxy/nodes, /galaxy/stage, /galaxy/window, /galaxy/near and /galaxy/territories, "
                            "fetched as the map's camera moves. Empty turns it off.", restart=True,
                            group="protection",
                            label="Map data limit",
                            summary="How fast one visitor may load map data while moving the map; "
                                    "empty turns the limit off.",
                            control="rate")
    health: str = opt("60 per minute", "Limit on /api/health. Empty turns it off.", restart=True,
                      group="protection",
                      label="Health check limit",
                      summary="How often one client may call the health check; empty turns the "
                              "limit off.",
                      control="rate")
    other: str = opt("300 per minute", "Limit on every other page, counted together. Empty turns it off.",
                     restart=True, group="protection",
                     label="Other pages limit",
                     summary="How often one visitor may load any other page, counted together; "
                             "empty turns the limit off.",
                     control="rate")


class RateLimit(_Section):
    default: str = opt("200 per day;50 per hour", "Flask-Limiter's default limit for the API.",
                       env="PLANETGEN_RATELIMIT_DEFAULT", restart=True, group="protection",
                       label="API request limit",
                       summary="How many API requests each client may make; empty turns the limit "
                               "off.",
                       control="rate")
    storage_uri: str = opt("", "Where limits count; empty means the Redis server in redis.url, memory:// counts "
                           "per process (right with exactly one worker).", env="PLANETGEN_RATELIMIT_STORAGE_URI",
                           restart=True, group="protection",
                           label="Where limits are counted",
                           summary="Empty counts in Redis; memory:// counts per process.",
                           advanced=True)
    pages: RatePages = RatePages()


class ApiLog(_Section):
    internal: Literal["all", "off"] = opt(
        "all", "Whether the calls the web pages make to the API in-process (about six for each page view) are written "
        "to the activity log as API lines: all, or off.", group="housekeeping",
        label="Log the pages' own API calls",
        summary="All API calls are logged by default; off leaves out the ones the web pages make to themselves.")


class ProxyFix(_Section):
    @field_validator("x_for", "x_proto", "x_host", mode="before")
    @classmethod
    def _not_a_bool(cls, value):
        if isinstance(value, bool):
            raise ValueError("must be a whole number of proxies (0 = off), not true or false")
        return value

    x_for: int = opt(0, "Reverse proxies to trust for X-Forwarded-For, the client address (0 = off).",
                     env="PLANETGEN_PROXY_FIX_X_FOR", restart=True, editable=False, group="proxy", ge=0,
                     unit="proxies", label="Proxies for the visitor's address",
                     summary="How many reverse proxies to trust for the visitor's address (0 means "
                             "none).")
    x_proto: int = opt(0, "Reverse proxies to trust for X-Forwarded-Proto, the scheme (0 = off).",
                       env="PLANETGEN_PROXY_FIX_X_PROTO", restart=True, editable=False, group="proxy", ge=0,
                       unit="proxies", label="Proxies for the scheme",
                       summary="How many reverse proxies to trust for http or https (0 means none).")
    x_host: int = opt(0, "Reverse proxies to trust for X-Forwarded-Host (0 = off).",
                      env="PLANETGEN_PROXY_FIX_X_HOST", restart=True, editable=False, group="proxy", ge=0,
                      unit="proxies", label="Proxies for the host name",
                      summary="How many reverse proxies to trust for the host name (0 means none).")


class PasswordHash(_Section):
    memory_kib: int = opt(19456, "Memory one admin password check uses, in KiB (Argon2id). 19456 (19 MiB) with "
                          "time_cost 2 and parallelism 1 is OWASP's minimum; raise it where the host has the memory "
                          "to spare. A stored hash made with other numbers is re-hashed at that admin's next sign-in.",
                          env="PLANETGEN_PASSWORD_HASH_MEMORY_KIB", unit="KiB", restart=True, editable=False,
                          group="security", ge=8,
                          label="Password hashing memory",
                          summary="How much memory checking one admin password takes; the default is "
                                  "OWASP's minimum.")
    time_cost: int = opt(2, "Passes over that memory for one admin password check (Argon2id).",
                         env="PLANETGEN_PASSWORD_HASH_TIME_COST", restart=True, editable=False, group="security",
                         ge=1,
                         unit="passes", label="Password hashing passes",
                         summary="How many passes one admin password check makes over its memory.")
    parallelism: int = opt(1, "Threads one admin password check uses (Argon2id).",
                           env="PLANETGEN_PASSWORD_HASH_PARALLELISM", restart=True, editable=False,
                           group="security", ge=1,
                           unit="threads", label="Password hashing threads",
                           summary="How many threads one admin password check uses.")


class Cache(_Section):
    max_mb: float = opt(1024, "How much disk the site's stored caches (today the Galaxy Map's nodes and stages) may use, in all. "
                        "When it is full the warm-up stops, the nodes it built ahead stay, and the nodes a visit "
                        "asked for are dropped, the least recently visited first (0 turns the disk caches off).",
                        env="PLANETGEN_CACHE_MAX_MB", unit="MB", group="speed", ge=0,
                        label="Disk cache size",
                        summary="How much disk the site's stored caches may use in all; 0 turns "
                                "them off.",
                        slider=(16, 1048576),
                        control="slider")
    warm_on_change: bool = opt(True, "Start the warm-up by itself after anything that changes what the site shows (a "
                               "generation run, the orbital update, an edit), so the caches are rebuilt without "
                               "waiting for the next update. It rebuilds only what the change touched.",
                               env="PLANETGEN_WARM_ON_CHANGE", group="speed",
                               label="Rebuild caches after changes",
                               summary="Start the cache warm-up by itself after anything that "
                                       "changes what the site shows.")
    warm_settle_seconds: float = opt(60, "How long the automatic warm-up waits after the first change before it "
                                     "starts, so a burst of edits is one warm-up and not one for each.",
                                     unit="seconds", group="speed", ge=0,
                                     label="Wait before rebuilding",
                                     summary="How long the warm-up waits after a change, so a burst "
                                             "of edits is one rebuild.",
                                     control="duration",
                                     advanced=True)

    _loose = field_validator("warm_on_change", mode="before")(_loose_bool)


class TileCache(_Section):
    dir: str = opt("", "Where the Galaxy Map's nodes are cached on disk; empty means /var/cache/planetgen/tiles.",
                   env="PLANETGEN_TILE_CACHE_DIR", editable=False, group="paths",
                   label="Map cache folder",
                   summary="Where the Galaxy Map's nodes are cached; empty uses "
                           "/var/cache/planetgen/tiles.")
    warm_zoom_floor: int = opt(9, "Finest block level of the zoom run `planetgen warm` builds from the opening view in to "
                               "the galaxy's centre, on top of the populated blocks it builds in full (9 is blocks of "
                               "nine shells, about 36 pc across; the run covers blocks of that size and bigger; "
                               "243 or more builds no run).",
                               unit="shells", group="speed", ge=1, le=177147,
                               label="Closest map block built ahead",
                               summary="How far in the map is prepared: blocks of this many shells; "
                                       "243 or more skips it.")
    prefetch_per_visit: int = opt(128, "Most nodes one request for nodes queues to be built in the background: "
                                  "the children of the nodes it asked for, then their parents and the parents' other "
                                  "children (0 turns the background build off).", unit="nodes", group="speed", ge=0,
                                  label="Nodes prefetched per visit",
                                  summary="Most nodes one visit queues to be built ahead in the "
                                          "background (0 turns it off).")
    prefetch_per_minute: int = opt(600, "Most nodes the background build queues in a minute, all requests together.",
                                   unit="nodes", group="speed", ge=0,
                                   label="Nodes prefetched per minute",
                                   summary="Most nodes the background build queues in a minute, all "
                                           "visitors together.")
    prefetch_max_queued: int = opt(8, "The background build queues nothing while this many of its jobs (16 nodes "
                                   "each) or other upkeep work already wait.", unit="jobs", group="speed", ge=1,
                                   label="Prefetch queue limit",
                                   summary="The background build waits while this many of its jobs "
                                           "are already queued.")


class PageCache(_Section):
    enabled: bool = opt(True, "Keep the API's public answers in memory, per WSGI process.",
                        env="PLANETGEN_PAGE_CACHE", restart=True, group="speed",
                        label="Page answer cache",
                        summary="Keep the API's public answers in memory so repeat requests are "
                                "fast.")
    max_entries: int = opt(2000, "Most answers kept.", restart=True, group="speed", ge=1,
                           unit="answers", label="Cached answers",
                           summary="The most answers kept in memory.",
                           advanced=True)
    max_mb: float = opt(64, "Most memory the answers may take.", restart=True, unit="MB", group="speed", gt=0,
                        label="Answer cache memory",
                        summary="The most memory the cached answers may use.",
                        advanced=True)
    stamp_seconds: float = opt(15, "How often the galaxy's content stamp is checked.", restart=True, unit="seconds",
                               group="speed", ge=0,
                               label="Change check interval",
                               summary="How often the galaxy's content stamp is checked.",
                               control="duration",
                               advanced=True)
    max_age_seconds: float = opt(300, "Nothing is kept longer than this.", restart=True, unit="seconds",
                                 group="speed", gt=0,
                                 label="Answer lifetime in memory",
                                 summary="No answer is kept in memory longer than this.",
                                 control="duration",
                                 advanced=True)
    shared_seconds: float = opt(21600, "How long an answer is kept in Redis, where every web process and the warm-up "
                                "share it. Each answer is served only while the data it read has not changed, so "
                                "this only bounds the memory.", restart=True, unit="seconds", group="speed", gt=0,
                                label="Answer lifetime in Redis",
                                summary="How long a shared answer is kept in Redis; it only bounds "
                                        "the memory.",
                                control="duration",
                                advanced=True)
    warm_max_objects: int = opt(5000, "Most sector pages, system pages and phenomenon pages the warm-up builds "
                                "(each, nearest the galaxy centre first); the rest are built on their first visit, "
                                "and the warm-up says how many it left out.", unit="objects", group="speed", ge=0,
                                label="Pages built ahead",
                                summary="Most sector, system and phenomenon pages the warm-up "
                                        "builds; the rest build on first visit.")
    stored_counts: bool = opt(True, "Table pages show stored totals and filter-menu counts, refreshed in the "
                              "background, instead of counting whole tables in the request (off: exact counts "
                              "at any price).", env="PLANETGEN_COUNT_CACHE", group="speed",
                              label="Use stored table counts",
                              summary="Table pages show stored totals, refreshed in the background, "
                                      "instead of counting on each request.")

    _loose = field_validator("enabled", "stored_counts", mode="before")(_loose_bool)


class Jobs(_Section):
    dir: str = opt("", "Where each background job's command lines, status and output are kept; empty means "
                   "/var/lib/planetGen/jobs.", env="PLANETGEN_JOBS_DIR", editable=False, group="paths",
                   label="Job files folder",
                   summary="Where each background job's command lines, status and output are kept.")
    keep: int = opt(20, "How many finished jobs are kept.", group="housekeeping", ge=1,
                    unit="jobs", label="Finished jobs kept",
                    summary="How many finished jobs stay in the job list.")
    python: str = opt("", "The interpreter that runs planetgen and planetgen.cli.reset for the jobs; empty means "
                      "the web app's own Python. Editing it from the web would be code execution, so it is file-only.",
                      env="PLANETGEN_PYTHON", editable=False, group="paths",
                      label="Python for jobs",
                      summary="The interpreter that runs the generation jobs; empty uses the web "
                              "app's own.")


class Parity(_Section):
    enabled: bool = opt(False, "Keep a parity file beside the database that can rebuild a damaged sector: every saved "
                        "or edited sector updates the parity of its group, which costs a few milliseconds a save and about "
                        "7 percent of the sectors' size in disk space. Off, saves write no parity and `planetgen "
                        "refresh-parity` builds it later.", env="PLANETGEN_PARITY", restart=True, group="housekeeping",
                        label="Sector repair file",
                        summary="Keep a parity file that can rebuild a damaged sector; saves get a little slower.")
    dir: str = opt("", "Where the parity file is kept; empty means /var/lib/planetGen/parity. Keep it off the disk the "
                   "database is on if you can: a parity file on the same failing disk protects nothing.",
                   env="PLANETGEN_PARITY_DIR", editable=False, group="paths",
                   label="Repair file folder",
                   summary="Where the sector repair file is kept; empty uses /var/lib/planetGen/parity.")
    group_size: int = opt(32, "Sectors in one parity group (G). Any `parity_blocks` damaged sectors of a group are rebuilt "
                          "from the rest, which are read to do it: a bigger group costs less disk and more reading. A "
                          "change applies to the groups made afterwards.", unit="sectors", restart=True, group="housekeeping",
                          ge=1, le=248, advanced=True, label="Sectors per repair group",
                          summary="How many sectors share one set of parity blocks.")
    parity_blocks: int = opt(2, "Parity blocks per group (m): how many damaged sectors of one group can be rebuilt. A "
                             "change applies to the groups made afterwards.", unit="blocks", restart=True,
                             group="housekeeping", ge=1, le=8, advanced=True, label="Repair blocks per group",
                             summary="How many damaged sectors of one group the repair file can rebuild.")


class Queue(_Section):
    interactive_workers: int = opt(
        1, "How many workers are kept for short jobs: edits, estimates and slow page answers. They "
        "are never counted in the bulk worker count, so a generation run on every core leaves them be, and they run "
        "at a higher priority than the bulk workers.", unit="workers", group="speed", ge=1, le=64,
        label="Workers for short jobs",
        summary="How many workers are reserved for short jobs; generation work never takes them.")
    maintenance_workers: int = opt(
        1, "How many workers run the cache warm-up, the node prefetch and the database check while a generation run "
        "is live. They are taken out of the generation worker count (a generation worker finishes its unit and "
        "leaves when they are needed), so the host's total does not grow; with no generation run live they may use "
        "the whole count.", unit="workers", group="speed", ge=1, le=64,
        label="Workers for upkeep during a run",
        summary="How many workers keep the caches and checks going while a generation run is live.")
    max_host_workers: int = opt(
        0, "The most worker processes planetgen may run on this host at once, every queue and command together: a "
        "start that would pass it is refused and said in the log. 0 sizes it from the host: 60 percent of its memory "
        "in workers of 300 MB each, no more than four per core, and never fewer than the cores plus four. It stops a "
        "burst of requests, or several commands starting workers together, from filling the memory and getting "
        "MySQL killed.", unit="workers", group="speed", ge=0, le=512,
        label="Most workers on this host",
        summary="The most worker processes planetgen runs at once on this host; 0 sizes it from memory and cores.")
    long_workers: int = opt(
        1, "How many workers run the jobs expected to take longer than the long job limit (a neighbourhood or a "
        "large settle, which can run for hours). They have a lane of their own so a long job can never hold the "
        "workers kept for short ones.", unit="workers", group="speed", ge=1, le=64,
        label="Workers for long jobs",
        summary="How many workers run long jobs, apart from the ones for short jobs.")
    long_job_seconds: float = opt(
        60, "A job expected to take longer than this goes to the long lane. A job with no estimate counts as short.",
        unit="seconds", group="speed", gt=0,
        label="Long job limit",
        summary="Jobs expected to take longer than this go to the long lane.",
        control="duration")
    unit_ceiling_seconds: float = opt(
        300, "The longest a unit of a generation stage is expected to run. A stage is cut into units no "
        "longer than this at the default galaxy size; a unit that runs longer is written to the run's log, naming "
        "the stage. A stage that is one unit by nature is not held to it. A unit this long is also the most a "
        "worker waits to hand its place to a lane that needs it.", unit="seconds", group="speed", gt=0,
        label="Unit length ceiling",
        summary="The longest one unit of a generation stage is meant to run; a longer one is logged.",
        control="duration")
    wait_minutes: float = opt(
        60, "How long a generation run, or an API job that generates sectors, waits in line when another generation "
        "run holds the machine, before it fails naming that run (0 does not wait).", unit="minutes", group="speed",
        ge=0,
        label="Longest wait for a run",
        summary="How long a run waits in line for the one in front of it before it gives up (0 means it "
                "does not wait).")
    max_waiting: int = opt(
        8, "How many entries the one waiting list holds: generation runs and API jobs that wait for the "
        "machine, and Generate page jobs that wait for the running job. One more than this is refused, naming the list.",
        unit="entries", group="speed", ge=1, le=100,
        label="Waiting list length",
        summary="How many runs and jobs can wait for their turn at once; one more is refused.")


class SettingsFile(_Section):
    backups_kept: int = opt(20, "How many backup archives are kept. One is made just before every settings save, "
                            "restore or import replaces the file, and by Back up now; the oldest is deleted once "
                            "there are more than this, but never the newest archive that holds every part.",
                            unit="backups", group="housekeeping", ge=1, le=200,
                            label="Settings backups kept",
                            summary="How many settings backup archives to keep; the oldest is "
                                    "deleted first.")


_HEX_COLOR = re.compile(r"^#[0-9a-fA-F]{6}$")
_ASSET_NAME = re.compile(r"^([0-9a-f]{16}\.(png|ico|jpg|svg))?$")
_TWITTER_HANDLE = re.compile(r"^(@[A-Za-z0-9_]{1,15})?$")
_LANGUAGE = re.compile(r"^[a-z]{2,3}(-[A-Za-z]{2,4})?$")

ROBOTS_DIRECTIVES = ("user-agent", "allow", "disallow", "sitemap", "crawl-delay", "host")
"""tuple[str]: The lines a robots.txt text may hold besides comments and blank lines (RFC 9309 and the common extras)."""


def robots_problems(text):
    """`[message]` for each line of a robots.txt text that is not a `Field: value` line with a known field, a comment
    or blank; empty when the text is fine. The first problems only, so one bad paste doesn't list a thousand."""
    problems = []
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        field, colon, value = line.partition(":")
        if not colon or field.strip().lower() not in ROBOTS_DIRECTIVES:
            problems.append(f"line {number} is not a robots.txt line (Field: value; known fields are "
                            f"{', '.join(ROBOTS_DIRECTIVES)})")
            if len(problems) >= 5:
                break
    return problems


class Site(_Section):
    tagline: str = opt("", "A short line shown under the home page's heading and used as the description of pages "
                       "that give none of their own.", group="site", max_length=120,
                       label="Tagline",
                       summary="A short line under the home page heading, and the description of "
                               "pages that have none.")
    language: str = opt("en", "The pages' language as a code (en, de, en-GB): the `lang` of every page and the Open "
                        "Graph locale.", group="site",
                        label="Language",
                        summary="The language code of the pages, such as en or en-GB.",
                        control="language")

    @field_validator("language")
    @classmethod
    def _language_code(cls, value):
        if not _LANGUAGE.match(value):
            raise ValueError("must be a language code such as en or en-GB")
        return value


class Web(_Section):
    title_template: str = opt("{title} - {site_name}", "How a page's title is written in the browser tab and in "
                              "search results; {title} is the page's own title and {site_name} the site name.",
                              group="site", max_length=120,
                              label="Page title format",
                              summary="How a page's title is written in the browser tab, using "
                                      "{title} and {site_name}.")
    theme_color_light: str = opt("#f6f7fb", "The browser's toolbar colour on phones for the light theme, as #rrggbb.",
                                 group="look",
                                 label="Browser colour, light theme",
                                 summary="The toolbar colour phones show with the light theme.",
                                 control="color")
    theme_color_dark: str = opt("#14151e", "The same for the dark theme.", group="look",
                                label="Browser colour, dark theme",
                                summary="The toolbar colour phones show with the dark theme.",
                                control="color")
    assets_dir: str = opt("", "Where uploaded icons and images are kept; empty means the assets folder in the state "
                          "folder.", env="PLANETGEN_ASSETS_DIR", editable=False, group="paths",
                          label="Uploaded images folder",
                          summary="Where uploaded icons and pictures are stored.")

    @field_validator("title_template")
    @classmethod
    def _title_has_title(cls, value):
        try:
            shown = value.format(title="X", site_name="Y")
        except (KeyError, IndexError, ValueError) as exc:
            raise ValueError(f"may only use {{title}} and {{site_name}} ({exc})") from exc
        if "X" not in shown:
            raise ValueError("must contain {title}")
        return value

    @field_validator("theme_color_light", "theme_color_dark")
    @classmethod
    def _color(cls, value):
        if not _HEX_COLOR.match(value):
            raise ValueError("must be a colour written #rrggbb")
        return value


class Icons(_Section):
    favicon_svg: str = opt("", "The site icon as an uploaded SVG; empty uses the built-in one. Set it by uploading on "
                           "the Images page.", group="look",
                           label="Site icon (SVG)",
                           summary="The icon browsers show in tabs; upload it on the Images page.",
                           control="asset")
    favicon_png: str = opt("", "A PNG or ICO icon for browsers that do not take an SVG; empty offers none. Set it by "
                           "uploading on the Images page.", group="look",
                           label="Site icon (PNG or ICO)",
                           summary="A fallback icon for browsers that do not read SVG; upload it on "
                                   "the Images page.",
                           control="asset")
    apple_touch_icon: str = opt("", "The 180 by 180 PNG iPhones and iPads use for a home-screen bookmark. Set it by "
                                "uploading on the Images page.", group="look",
                                label="Home-screen icon",
                                summary="The icon iPhones and iPads use when the site is added to "
                                        "the home screen.",
                                control="asset")
    manifest_enabled: bool = opt(False, "Serve /site.webmanifest, which names the site, its icon and its theme colour "
                                 "for browsers that can install it.", group="look",
                                 label="Installable site",
                                 summary="Let browsers offer to install the site as an app.")

    _loose = field_validator("manifest_enabled", mode="before")(_loose_bool)

    @field_validator("favicon_svg", "favicon_png", "apple_touch_icon")
    @classmethod
    def _asset(cls, value):
        if not _ASSET_NAME.match(value):
            raise ValueError("is the name of an uploaded file; upload it on the Images page")
        return value


class Seo(_Section):
    indexing: Literal["off", "top_pages", "all"] = opt(
        "top_pages", "Whether search engines may index the site: off (robots.txt disallows everything and every page "
        "says noindex), top_pages (the home, class, species, polity, galaxy and list pages), or all (also the "
        "system, sector and phenomenon pages).", group="search",
        choice_labels={"off": "Off", "top_pages": "Main pages only", "all": "All pages"},
        label="Search engine indexing",
        summary="Which pages search engines may list: none, the main pages, or all of them.")
    detail_pages: Literal["noindex", "index"] = opt(
        "noindex", "Whether the generated system, sector and phenomenon pages may be indexed when indexing is "
        "top_pages. Millions of generated pages look like scaled content to search engines, so they stay out by "
        "default; a small galaxy can switch this to index.", group="search",
        choice_labels={"noindex": "Keep out of search", "index": "Allow in search"},
        label="Index generated pages",
        summary="Whether the generated system, sector and phenomenon pages may be indexed.")
    description: str = opt("", "The description search engines show for pages that give none of their own (about 155 "
                           "characters are shown); empty uses the tagline.", group="search", max_length=300,
                           label="Search description",
                           summary="The snippet search engines show for pages without their own; "
                                   "empty uses the tagline.")
    keywords: List[str] = opt([], "Keywords for the keywords meta tag. Google has ignored it since 2009; some "
                              "smaller engines still read it.", group="search",
                              label="Keywords",
                              summary="Words for the keywords tag; most search engines ignore them.")
    json_ld: bool = opt(True, "Add a WebSite block (schema.org) to the home page.", group="search",
                        label="Structured data",
                        summary="Add a website block to the home page so search engines understand "
                                "the site.")
    robots_extra: str = opt("", "Lines added to the end of the generated robots.txt.", group="search", multiline=True,
                            max_length=64 * 1024,
                            label="Extra robots.txt lines",
                            summary="Lines added to the end of the generated robots.txt.",
                            advanced=True)
    robots_override: str = opt("", "Replaces the generated robots.txt entirely when not empty. A mistake here can "
                               "hide the whole site from search engines.", group="search", multiline=True,
                               max_length=500 * 1024,
                               label="Replace robots.txt",
                               summary="Replaces the whole generated robots.txt; a mistake can hide "
                                       "the site from search engines.",
                               advanced=True)
    block_ai_crawlers: bool = opt(False, "Disallow the crawlers that collect text to train AI models (GPTBot, "
                                  "ClaudeBot, CCBot, Google-Extended, Applebot-Extended) in robots.txt.",
                                  group="search",
                                  label="Block AI crawlers",
                                  summary="Ask the crawlers that collect text to train AI models to "
                                          "stay away.")
    sitemap_enabled: bool = opt(True, "Serve /sitemap.xml and list it in robots.txt.", group="search",
                                label="Sitemap",
                                summary="Publish /sitemap.xml and list it in robots.txt.")
    sitemap_max_urls: int = opt(50000, "Most addresses in one sitemap file (the limit is 50,000).", unit="addresses",
                                group="search", ge=1, le=50000,
                                label="Sitemap size",
                                summary="The most addresses in one sitemap file (up to 50,000).",
                                advanced=True)
    sitemap_cache_hours: int = opt(24, "How long a built sitemap is reused.", unit="hours", group="search", ge=1,
                                   le=168,
                                   label="Sitemap refresh",
                                   summary="How long a built sitemap is reused before it is built "
                                           "again.",
                                   control="duration",
                                   advanced=True)
    security_txt: str = opt("", "The text of /.well-known/security.txt (RFC 9116 needs a Contact and an Expires "
                            "line); empty serves none.", group="search", multiline=True, max_length=16 * 1024,
                            label="security.txt",
                            summary="The text served at /.well-known/security.txt; empty serves "
                                    "none.",
                            advanced=True)
    humans_txt: str = opt("", "The text of /humans.txt; empty serves none.", group="search", multiline=True,
                          max_length=16 * 1024,
                          label="humans.txt",
                          summary="The text served at /humans.txt; empty serves none.",
                          advanced=True)
    privacy_note: str = opt("", "A short note shown on the account page about what the site keeps.", group="search",
                            multiline=True, max_length=4 * 1024,
                            label="Privacy note",
                            summary="A short note on the account page about what the site keeps.")

    _loose = field_validator("json_ld", "block_ai_crawlers", "sitemap_enabled", mode="before")(_loose_bool)

    @field_validator("robots_extra", "robots_override")
    @classmethod
    def _robots(cls, value):
        problems = robots_problems(value)
        if problems:
            raise ValueError("; ".join(problems))
        return value

    @field_validator("keywords")
    @classmethod
    def _keywords(cls, value):
        cleaned = [word.strip() for word in value if word.strip()]
        if len(", ".join(cleaned)) > 300:
            raise ValueError("the keywords together may be at most 300 characters")
        return cleaned


class OpenGraph(_Section):
    enabled: bool = opt(True, "Add Open Graph and Twitter tags to the public pages, so a link shows a title, a "
                        "description and a picture when shared.", group="previews",
                        label="Link previews",
                        summary="Show a title, description and picture when a link to the site is "
                                "shared.")
    site_name: str = opt("", "The site name in link previews; empty uses the site name.", group="previews",
                         max_length=100,
                         label="Name in previews",
                         summary="The site name shown in link previews; empty uses the site name.")
    image_alt: str = opt("", "A description of the default preview picture for people who cannot see it.",
                         group="previews", max_length=200,
                         label="Preview picture description",
                         summary="Describes the preview picture for people who cannot see it.")
    default_image: str = opt("", "The picture link previews show (a PNG or JPEG about 1200 by 630); empty uses the "
                             "built-in one. Set it by uploading on the Images page.", group="previews",
                             label="Preview picture",
                             summary="The picture link previews show; upload it on the Images page.",
                             control="asset")
    twitter_card: Literal["summary", "summary_large_image"] = opt(
        "summary_large_image", "How X (Twitter) draws a preview: a large picture or a small one.", group="previews",
        choice_labels={"summary": "Small picture", "summary_large_image": "Large picture"},
        label="X (Twitter) card style",
        summary="Whether X draws a large or a small preview picture.")
    twitter_site: str = opt("", "The site's X (Twitter) account, as @name; empty names none.", group="previews",
                            label="X (Twitter) account",
                            summary="The site's X account as @name; empty names none.")
    generated_cards: bool = opt(False, "Draw a preview picture per system and sector instead of one for the whole "
                                "site. Off until a card drawer is installed; every page then uses the "
                                "default picture.", group="previews",
                                label="Preview picture per page",
                                summary="Draw a preview picture for each system and sector; off "
                                        "until a card drawer is installed.")
    card_cache_mb: float = opt(100, "How much disk the drawn preview pictures may use.", unit="MB", group="previews",
                               ge=0,
                               label="Preview picture disk use",
                               summary="How much disk the drawn preview pictures may use.",
                               slider=(10, 10240),
                               control="slider")
    card_cache_dir: str = opt("", "Where drawn preview pictures are cached; empty means the state folder.",
                              editable=False, group="paths",
                              label="Preview picture folder",
                              summary="Where drawn preview pictures are kept.")

    _loose = field_validator("enabled", "generated_cards", mode="before")(_loose_bool)

    @field_validator("default_image")
    @classmethod
    def _image(cls, value):
        if not _ASSET_NAME.match(value):
            raise ValueError("is the name of an uploaded file; upload it on the Images page")
        return value

    @field_validator("twitter_site")
    @classmethod
    def _handle(cls, value):
        if not _TWITTER_HANDLE.match(value):
            raise ValueError("must be an @name of up to 15 letters, digits or underscores")
        return value


class Settings(_Section):
    """Every `config.json` option."""

    config_version: int = opt(CONFIG_VERSION, "The settings file format version; a file without one is version 1.",
                              editable=False, group="file", ge=1,
                              unit="format version", label="Settings file version",
                              summary="Which format version this settings file is in; the program "
                                      "keeps it up to date.")
    site_name: str = opt("planetGen", "Display name for this deployment, shown in the header, the home page "
                         "heading, every page title and the footer.", group="site", min_length=1,
                         label="Site name",
                         summary="The name shown in the header, the home page, page titles and the "
                                 "footer.")
    base_url: str = opt("http://localhost/", "The address this site is reached at, as an absolute http or https "
                        "address with a trailing slash. Canonical links, link previews, the sitemap and emailed links "
                        "are built from it, never from the request's Host header.", role="owner", group="site",
                        label="Site address",
                        summary="The public address of the site, used for canonical links, previews "
                                "and the sitemap.")

    @field_validator("base_url")
    @classmethod
    def _base_url(cls, value):
        parts = urlsplit(value)
        if parts.scheme not in ("http", "https") or not parts.netloc or parts.query or parts.fragment \
                or parts.username or parts.password:
            raise ValueError("must be an absolute http or https address such as https://example.org/ (no query, "
                             "fragment or password)")
        return value if value.endswith("/") else value + "/"

    api_base_url: str = opt("", "Base URL of the Flask API's /api mount point (for example https://example.org/api). "
                            "The update asks it whether the site answers, and web/lib/apiclient.py calls it outside "
                            "the Flask app. There is no default: nothing guesses a loopback address, so it must be "
                            "set for either to work.", env="PLANETGEN_API_BASE_URL",
                            restart=True, editable=False, group="paths",
                            label="API address",
                            summary="Where the API lives; the update's health check asks it. No default.")

    @field_validator("api_base_url")
    @classmethod
    def _api_base_url(cls, value):
        if not value:
            return value
        parts = urlsplit(value)
        if parts.scheme not in ("http", "https") or not parts.netloc or parts.query or parts.fragment \
                or parts.username or parts.password:
            raise ValueError("must be empty or an absolute http or https address such as https://example.org/api "
                             "(no query, fragment or password)")
        return value.rstrip("/")
    debug: bool = opt(False, "Debug mode: every decision the generator makes, every "
                      "random roll, every SQL statement, queue event and request, written to the console of each "
                      "planetgen command (stdout) and to log_file. The output is large.", env="PLANETGEN_DEBUG",
                      restart=True, editable=False, group="logging",
                      label="Debug log",
                      summary="Write a very detailed log of everything the generator does; it grows "
                              "fast.")
    debug_areas: str = opt("", "Limit debug mode to some areas of output: decision, draw, sql, queue, cache, web, lib "
                           "(comma separated; a leading - leaves one out). Empty writes everything but branch, "
                           "which logs every conditional jump of the program's own code and is written only "
                           "when named.",
                           env="PLANETGEN_DEBUG_AREAS", restart=True, editable=False, group="logging",
                           label="Debug areas",
                           summary="Which kinds of detail debug mode writes; empty is all of them.")
    log_file: str = opt(logpaths.DEFAULT_LOG_FILE, "Where the debug log goes.",
                        env="PLANETGEN_LOG_FILE", restart=True, editable=False, group="logging",
                        label="Debug log file",
                        summary="The file the debug log is written to.")
    record_file: str = opt("", "A file every log line of a planetgen command is also written to as a record: one per "
                           "line, the time in UTC. Empty writes none; the --record-file option does the same for "
                           "one command.",
                           env="PLANETGEN_RECORD_FILE", restart=True, editable=False, group="logging",
                           label="Record file",
                           summary="A file that gets every log line as a machine-readable record, for a tool to read.")
    record_format: Literal["logfmt", "json"] = opt(
        "logfmt", "The form of the record file's lines: logfmt (key=value) or json (one object per line).",
        env="PLANETGEN_RECORD_FORMAT", restart=True, editable=False, group="logging",
        choice_labels={"logfmt": "Key and value pairs", "json": "JSON lines"},
        label="Record file form",
        summary="Whether the record file is written as key=value lines or as JSON lines.")
    log_dir: str = opt("", "The folder of the always-on activity log; empty means /var/log/planetgen on Linux, "
                       "/Library/Logs/planetgen on macOS.",
                       env="PLANETGEN_LOG_DIR", restart=True, editable=False, group="logging",
                       label="Activity log folder",
                       summary="The folder that holds the always-on activity log.")
    log_rotation: Literal["auto", "system", "app"] = opt(
        logpaths.DEFAULT_LOG_ROTATION, "How the activity log is rotated: system (logrotate or newsyslog), app (the "
        "program, 100 MB per file, 30 old copies) or auto (system where an installed rotation config exists, else "
        "app).", restart=True, editable=False, group="logging",
        choice_labels={"auto": "Automatic", "system": "The system (logrotate)", "app": "The program"},
        label="Log rotation",
        summary="Who trims old activity logs: the system, the program, or whichever is installed.")
    settings: SettingsFile = SettingsFile()
    site: Site = Site()
    web: Web = Web()
    icons: Icons = Icons()
    seo: Seo = Seo()
    og: OpenGraph = OpenGraph()
    mysql: MySQL = MySQL()
    control_database: str = opt("planetgen_control", "The MySQL schema holding admin identities, sessions, API keys "
                                "and the audit log. Never listed or selectable as a galaxy.",
                                env="PLANETGEN_CONTROL_DATABASE", restart=True, editable=False, group="database",
                                label="Control database",
                                summary="The database that holds admin accounts, sessions, API keys "
                                        "and the audit log.")
    redis: Redis = Redis()
    queue: Queue = Queue()
    ratelimit: RateLimit = RateLimit()
    api_log: ApiLog = ApiLog()
    login_allowlist: List[str] = opt(
        [], "Client addresses or networks (CIDR) never locked out after failed logins; loopback never is either. "
        "Environment: comma- or space-separated.", env="PLANETGEN_LOGIN_ALLOWLIST", restart=True, role="owner",
        group="security",
        label="Never-locked-out addresses",
        summary="Client addresses or networks that are never locked out after failed logins.")
    admin_cookie_insecure: bool = opt(False, "Send the admin session cookie over plain HTTP. Only for local "
                                      "development without TLS; a production deployment must never set it.",
                                      env="PLANETGEN_ADMIN_COOKIE_INSECURE", restart=True, editable=False,
                                      group="security",
                                      label="Admin login over plain HTTP",
                                      summary="Only for local development without TLS; never turn "
                                              "it on in production.")
    secret_key: str = opt("", "Signs the CSRF tokens on the pages' forms. Empty makes one and keeps it in the state folder "
                          "(mode 0600), so every worker and restart shares it.", env="PLANETGEN_SECRET_KEY", restart=True,
                          secret=True, editable=False, group="security",
                          label="Form signing key",
                          summary="Signs the security tokens on forms; empty makes a random one at "
                                  "each start.")
    password_hash: PasswordHash = PasswordHash()
    proxy_fix: ProxyFix = ProxyFix()
    cache: Cache = Cache()
    tile_cache: TileCache = TileCache()
    page_cache: PageCache = PageCache()
    jobs: Jobs = Jobs()
    parity: Parity = Parity()

    _loose = field_validator("debug", "admin_cookie_insecure", mode="before")(_loose_bool)

    @field_validator("login_allowlist", mode="before")
    @classmethod
    def _split_allowlist(cls, value):
        if isinstance(value, str):
            return value.replace(",", " ").split()
        return value


class SettingsError(ValueError):
    """A settings file or environment variable that doesn't fit the model;
    the message names each field."""


# ---------------------------------------------------------------------
# Walking the model
# ---------------------------------------------------------------------

def _model_of(annotation):
    if isinstance(annotation, type) and issubclass(annotation, BaseModel):
        return annotation
    return None


def iter_fields(model=Settings, prefix=""):
    """
    Yields `(dotted path, FieldInfo, metadata dict)` for every leaf option
    of `model` (sections are walked into). The metadata always has every
    `x-` key (`opt`'s defaults), `description` and the annotation.
    """
    for name, info in model.model_fields.items():
        path = f"{prefix}{name}"
        section = _model_of(info.annotation)
        if section is not None:
            yield from iter_fields(section, f"{path}.")
            continue
        meta = dict(info.json_schema_extra or {})
        meta["description"] = info.description or ""
        meta["annotation"] = info.annotation
        yield path, info, meta


def field_default(info):
    return info.get_default(call_default_factory=True)


def type_name(annotation):
    """`"integer"`, `"text"`, `"list of text"`, ... for the docs table."""
    origin = get_origin(annotation)
    if origin is Literal:
        return " or ".join(f'"{a}"' for a in get_args(annotation))
    if origin in (list, List):
        return "list of text"
    return {bool: "true or false", int: "whole number", float: "number", str: "text"}.get(annotation, str(annotation))


def env_names():
    """Every environment variable the model names, as `{variable: path}`."""
    return {meta["x-env"]: path for path, _info, meta in iter_fields() if meta.get("x-env")}


def env_name(path):
    """The environment variable the model names for option `path`."""
    for found, _info, meta in iter_fields():
        if found == path:
            return meta["x-env"]
    raise KeyError(path)


# ---------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------

DEFAULT_STATE_DIR = "/var/lib/planetGen/state"
"""str: The web-owned state folder (the installers give it to Apache's user; `planetgen.admin.maintenance`
keeps its switch there too). The install folder above it is root's, so a file beside `config.json` could be
read by the web app but never replaced."""


def state_dir():
    """The web-owned state folder: `PLANETGEN_STATE_DIR`, else `/var/lib/planetGen/state`."""
    return os.environ.get("PLANETGEN_STATE_DIR") or DEFAULT_STATE_DIR


def _read_json_object(path, label):
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except FileNotFoundError:
        return {}
    except (OSError, ValueError) as exc:
        raise SettingsError(f"{label} ({path}) can't be read: {exc}") from exc
    if not isinstance(data, dict):
        raise SettingsError(f"{label} ({path}) must contain a JSON object, not {type(data).__name__}")
    return data


def migrate(data):
    """
    Brings a file's dict up to `CONFIG_VERSION` (no migrations exist yet:
    version 1 is the first).

    Raises:
        SettingsError: The file is from a newer release.
    """
    version = data.get("config_version", 1)
    if not isinstance(version, int) or isinstance(version, bool) or version < 1:
        raise SettingsError(f"config_version must be a whole number 1 or above, not {version!r}")
    if version > CONFIG_VERSION:
        raise SettingsError(f"config_version {version} is newer than this release understands ({CONFIG_VERSION}); "
                            "update planetGen")
    return data


def _merge(base, overrides):
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _merge(base[key], value)
        else:
            base[key] = copy.deepcopy(value)  # never share a layer's nested dict with the merge


def _set_path(data, path, value):
    parts = path.split(".")
    for part in parts[:-1]:
        if not isinstance(data.get(part), dict):
            data[part] = {}
        data = data[part]
    data[parts[-1]] = value


def _environment_layer(environ):
    layer = {}
    for path, _info, meta in iter_fields():
        name = meta.get("x-env")
        if not name or name not in environ:
            continue
        value = environ[name]
        if value == "" and not meta.get("x-env-empty"):
            continue
        _set_path(layer, path, value)
    return layer


def build(config=None, environ=None, strict=False):
    """
    Validates the layers into a `Settings`: defaults, `config` (the
    `config.json` dict), then the `PLANETGEN_*` variables in `environ`.

    Args:
        strict (bool): Unknown keys in `config` are errors (`config check`)
            instead of ignored.

    Raises:
        SettingsError: Naming every field that doesn't fit.
    """
    environ = os.environ if environ is None else environ
    merged = {}
    problems = []
    if config is not None:
        try:
            layer = migrate(dict(config))
        except SettingsError as exc:
            raise SettingsError(f"config.json: {exc}") from exc
        if strict:
            known = {path for path, _i, _m in iter_fields()}
            problems.extend(f"config.json: {message}" for message in _unknown_paths(layer, known))
        _merge(merged, layer)
    _merge(merged, _environment_layer(environ))
    if problems:
        raise SettingsError("\n".join(problems))
    try:
        return Settings.model_validate(merged)
    except ValidationError as exc:
        lines = [f"{path}: {message}" for path, message in _error_pairs(exc)]
        raise SettingsError("; ".join(lines) if len(lines) < 4 else "\n".join(lines)) from exc


def _error_pairs(exc):
    """`[(dotted path, message)]` of a `ValidationError`, the path cut back to the option it names (a list
    option's element error is the list's)."""
    leaves = {path for path, _i, _m in iter_fields()}
    pairs = []
    for error in exc.errors():
        parts = [str(part) for part in error["loc"]]
        path = next((".".join(parts[:n]) for n in range(len(parts), 0, -1) if ".".join(parts[:n]) in leaves),
                    ".".join(parts))
        pairs.append((path, error["msg"]))
    return pairs


def field_errors(config=None, environ=None):
    """
    Validates `config` like `build` but without raising and without the unknown-key check: `(settings, [])` when
    it fits, else `(None, [(dotted option path, message), ...])`. The environment is left out unless `environ`
    is given. The writer uses it to put each message at its option.
    """
    environ = {} if environ is None else environ
    merged = {}
    if config:
        _merge(merged, migrate(dict(config)))
    _merge(merged, _environment_layer(environ))
    try:
        return Settings.model_validate(merged), []
    except ValidationError as exc:
        return None, _error_pairs(exc)


def _unknown_paths(layer, known, prefix=""):
    sections = {path.rsplit(".", 1)[0] for path in known if "." in path}
    out = []
    for key, value in layer.items():
        path = f"{prefix}{key}"
        if path in known:
            continue
        if isinstance(value, dict) and (path in sections or any(k.startswith(f"{path}.") for k in known)):
            out.extend(_unknown_paths(value, known, f"{path}."))
        else:
            out.append(f"{path} is not a setting")
    return out


_lock = threading.Lock()
_cache = {"key": None, "settings": None}


def _stamp(path):
    try:
        stat = os.stat(path)
    except OSError:
        return None
    return (stat.st_mtime_ns, stat.st_size)


def _environment_stamp(environ):
    return tuple(sorted((name, environ[name]) for name in env_names() if name in environ))


def get_settings(environ=None):
    """
    The validated settings. Re-read only when `config.json` or a variable
    the model names changes (a stat costs about a microsecond), so every
    option that is read per use stays live.

    Raises:
        SettingsError: A file or variable doesn't fit the model.
    """
    environ = os.environ if environ is None else environ
    config_path = logpaths.CONFIG_PATH
    key = (config_path, _stamp(config_path), _environment_stamp(environ))
    with _lock:
        if _cache["key"] == key:
            return _cache["settings"]
        settings = build(_read_json_object(config_path, "config.json"), environ)
        _cache["key"], _cache["settings"] = key, settings
        return settings


def reset_cache():
    """Forgets the cached settings (tests)."""
    with _lock:
        _cache["key"] = _cache["settings"] = None


# ---------------------------------------------------------------------
# What the model writes: the docs table, the example file, the schema
# ---------------------------------------------------------------------

def example_dict():
    """`config.json.example`'s content: every default of the options
    `config.json` may hold, in the model's order."""
    return Settings().model_dump(mode="json")


def example_json():
    return json.dumps(example_dict(), indent=2) + "\n"


def schema_json():
    """`config.schema.json`: the model's JSON Schema."""
    return json.dumps(Settings.model_json_schema(), indent=2, sort_keys=True) + "\n"


def _format_default(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    if value == "":
        return "empty"
    if isinstance(value, list):
        return "empty list" if not value else ", ".join(str(v) for v in value)
    return str(value)


def docs_table():
    """The generated option table for `docs/reference/config.md` (Markdown)."""
    rows = ["| Option | Setting | Type | Default | Applies | Edited | Environment | What it does |",
            "|---|---|---|---|---|---|---|---|"]
    roles = {"admin": "any admin", "owner": "the Owner"}
    for path, info, meta in iter_fields():
        applies = "after a restart" if meta["x-restart"] else "at once"
        edited = roles[meta["x-min_role"]] if meta["x-editable"] else "command line only"
        unit = f" ({meta['x-unit']})" if meta.get("x-unit") else ""
        secret = " **Secret.**" if meta["x-secret"] else ""
        env = f"`{meta['x-env']}`" if meta.get("x-env") else "none"
        rows.append(f"| `{path}` | {meta['x-label']} | {type_name(meta['annotation'])}{unit} | {_format_default(field_default(info))} | "
                    f"{applies} | {edited} | {env} | {meta['description']}{secret} |")
    return "\n".join(rows) + "\n"


def _lookup(settings, path):
    node: Any = settings
    for part in path.split("."):
        node = getattr(node, part)
    return node
