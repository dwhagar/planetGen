# API design standards: versions, key scopes, uploads, logging, recipes

The rules the API needs before remote generation (API.3) and generate-by-recipe (API.18) are built: how the API is versioned and how a client learns it is too old, what an API key may do, how large and how fast an upload may be, how an upload is staged and checked, how every call is logged, and how a recipe is validated. It records what the repo does today, the measurements behind the numbers, and a recommended design for each item. Nothing here is built.

Informs: API.3, API.4, API.5, API.6, API.7, API.8, API.9, API.10, API.11, API.12, API.13, API.14, API.15, API.16, API.17, API.18, API.19, ADM.13

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [S] seen in a source (URL in section 11), [C] computed or measured in the research environment, [R] recalled and unconfirmed (listed in section 11). The research environment could only read search-result text and raw files from GitHub and PyPI, not papers or RFCs, so every outside-standards claim is [S] from raw IETF, Apache, tus, Pydantic or oasdiff source text, or [R]. Absolute MB/s and millisecond figures came from a slow sandbox CPU and MariaDB 10.11 (CI targets 11.4); treat them as relative.

Related: [reproducible-galaxies.md](reproducible-galaxies.md) (seed, version key, API.16 and API.17), [login-brute-force-protection.md](login-brute-force-protection.md) (lockouts, two-factor), [work-queue-audit.md](work-queue-audit.md) (which API paths run on the queue), [library-migration.md](library-migration.md), [object-ids.md](object-ids.md), [../api.md](../api.md) (the API as built).

## 1. Decisions already taken

- Remote generation (Boss, 2026-10-01 19:22Z to 19:32Z, API.3): the local machine generates in memory and needs no database; name indexes and id state are downloaded once; the server reserves id blocks and claims sectors for the run; uploads and downloads are compressed ("gzip, bz2 or similar"); the run is fully resumable (the client caches every call on disk until the server confirms it; the server buffers, stages and writes only complete units, a star system or a sector); the local API version must match the server's (API.4, API.5); only an admin's key can upload, user-level keys can read but never upload (API.6); upload limits are planned separately (API.7); reserved ids and claimed sectors stay reserved until the upload finishes or an admin clears it (ADM.13); the server checks every upload and has the final say (API.8).
- API.12 (Plan 2026-10-07, superseded 2026-10-08): the download carries the star and sector name registries and the naming key.
- Keys (Boss, 2026-10-02 01:31Z, API.6): every key is attached to the account that created it, and requests made with it act as that account.
- Logging (Boss, 2026-10-02 01:31Z, API.15): by default every API call is logged with the user (key owner, signed-in user, or `god` for the console), how it came in (API key, web session or console) and the HTTP response code; the API docs say HTTP codes are used for parity with web server logs.
- Recipes (Boss, 2026-10-07 11:47Z, API.18, API.19): one adaptable JSON system for sectors, systems, planets, moons, phenomena and a galaxy, any field fixed, ranged or random; validation errors are JSON with details and the command's log output; 400 for a gibberish request, 422 for a validation failure; a galaxy is built region by region.
- The queue (Boss, 2026-10-01 22:13Z, plan approved 2026-10-07): everything that communicates with the API uses the work queue whenever possible; long work answers `202` with a job id, short generation is queued and waited on, plain row writes stay in the request ([work-queue-audit.md](work-queue-audit.md); built for the routes in [../api.md](../api.md#queued-jobs)).
- Libraries (Boss, 2026-10-07): request bodies are Pydantic models (ADM.21, built); request limits and login lockouts run on Flask-Limiter with Redis (SEC.30, done, PR #643); two-factor sign-in runs on pyotp and segno (SEC.29, done, PR #612). Both SEC items are verified finished (SEC rows of [todo-number-map.md](todo-number-map.md); CHANGELOG.md entries for the Redis rate limits and the segno QR codes) and appear below only as constraints.

## 2. Summary of the recommendation

Version the remote-client contract only, by hand against a CI rule (4). Upload the generator's object graph in a versioned envelope, not table rows (6.1). Four key scopes in a join table, old keys grandfathered as `admin` (5). Gzip batches of 1 to 2 MiB with a 4 MiB server limit, per-path proxy exceptions and a bounded decompressor (6). Idempotent numbered batches staged on a disk spool with bookkeeping tables, not mirror tables (6.5, 7.3). Log every call to the activity log first (8). Recipes as Pydantic discriminated unions, 400 for gibberish and 422 for validation failures (9). Two bugs to fix on the way: deep JSON nesting gives a 500 (9.5), and Flask-Limiter sends `Retry-After` on 200 responses (6.7).

## 3. What the repo does today (checked 2026-10-09)

| Fact | Where | Consequence |
|---|---|---|
| Keys are `pg_` + `secrets.token_urlsafe(32)` (256 bits), stored as SHA-256 hex, `UNIQUE(key_hash)` | `admin/auth.py`, `control_schema.sql` | Storage is right; scopes, expiry and a prefix are missing. |
| `validate_api_key` runs `UPDATE ... last_used_at` and `commit()` on every authenticated request | `admin/auth.py` | A write per call; throttle it. |
| `admin_api_keys`: id, admin_user_id, label, key_hash, created_at, last_used_at, revoked_at | `control_schema.sql` | `CONTROL_SCHEMA_VERSION = 10` (`db/store.py`; v10 is DB.13), so the scope migration is control v11. API.9's text says "v8", which is out of date. Control migrations are still legacy steps in `store.py`; only the galaxy schema (v65) moved to Alembic revisions after v61. |
| Every write route is `require_admin`, key or session; `session_only` for key, credential, 2FA and logout routes | `api/authz.py` | A scope check slots in beside `fresh` and `session_only`. |
| Read routes are public, 200/day and 50/hour per IP; writes 10/minute | `api/limiter.py`, `api/routes.py` | A key-holding client needs key-based limits or the 50/hour default ruins it. |
| Bodies are Pydantic models; failure is `400` with `errors:[{field,message}]` | `api/schemas.py` | Recipes get 422; old routes stay 400 (9.3). |
| `MAX_CONTENT_LENGTH = 2 MiB`, enforced early by `_reject_oversized_body` | `api/config.py`, `web/app.py` | A per-route override works (6.3). |
| Proxy body limits: nginx `1m`; Caddy `max_size 1MB` | `examples/nginx`, `macos`, `caddy` | Below the app's 2 MiB; upload paths need an exception in every file. |
| One WSGI process, 5 threads, 60 s request timeout | `examples/apache/planetgen.conf.example` | The concurrency budget is tiny. |
| Activity log categories `AUTH`, `AUTHZ`, `DB`, `GEN`; other categories are dropped | `admin/activity_log.py` | API.15 must add one. |
| `StarSystem.to_dict()` and `from_dict()` exist and carry `schema_version` | `generation/system.py` | The natural upload unit. |
| Python 3.9 gets Flask-Limiter 3.11.0 and limits 4.2; 3.10+ gets 4.1.1 and 5.8.0 | `requirements.lock` | `cost` callables and `key_func` exist in both; confirm on 3.11.0 in the 3.9 CI job. |

## 4. Versioning and compatibility (API.3, API.4, API.5, API.16, API.17)

### 4.1 Who can be skewed, and which scheme

Three callers: the browser pages (same release, in-process through `web/transport.py`), `generate.py --remote` (installed elsewhere, can lag or lead the server), and third-party scripts. Only the second needs a handshake; third-party reads stay additive-only under a documented promise.

A `/api/v1/` tree duplicates 101 routes for one skewed client. Date-based, account-pinned versions (Stripe `Stripe-Version`, GitHub `X-GitHub-Api-Version: 2022-11-28`) need per-account pinning and a converter chain built for thousands of integrators [R]. Content negotiation is awkward in `curl` and Flask for no benefit. A custom header carrying `MAJOR.MINOR` fits (no `X-` prefix, which RFC 6648 deprecated [R]); PATCH has no meaning for a contract.

### 4.2 Handshake (API.5)

- The client sends `PlanetGen-API-Version: 1.4` on every request; the server sends its own on every response.
- `GET /api/version` (public, DB-free, a looser limit than `/api/health`'s 60 per minute) returns `{"api_version":"1.4","min_client":"1.2","max_client":"1.4","release":"7.380.0","version_key":"<22 hex>","galaxy_schema":65,"remote":{"upload":true}}`. `/api/health` also gets `api_version`, but the check must not depend on it (health is rate-limited and hits the database).
- Acceptance: same MAJOR and client MINOR <= server MINOR is accepted as is. A lower MAJOR is accepted only if a converter exists and its window is open. This answers API.5's "how many older versions": the current major, plus the previous major for one release cycle after a major bump, for upload payloads only. Converters exist only for what clients write, never for reads.
- Refusal: `400 {"error":"api_version_unsupported","your_version":"0.9","min":"1.2","max":"1.4","install":"update to release 7.380.x or newer"}`; `410 Gone` past a sunset. `426` is avoided: in RFC 9110 it signals a protocol switch and needs an `Upgrade` header [R].
- `generate.py --remote` calls `/api/version` before reading config or starting a worker and prints both versions and the fix.
- Deprecation: `Sunset` (RFC 8594) and `Deprecation` (RFC 9745, March 2025, a structured-field date such as `@1688169599`) with `Link` relations `deprecation` and `sunset` [R]. Send them only to clients whose declared major is older than current, and log the client version (API.15) to see who still uses an old major.

### 4.3 Machine-readable record: OpenAPI generation (PyPI facts, 2026-10-09)

Pydantic is at 2.14.0 (2026-10-08, Python >=3.10; the last for 3.9 is 2.13.5, which the lock pins) and `model_json_schema()` gives JSON Schema 2020-12, which OpenAPI 3.1 uses [R]; discriminated unions emit the OpenAPI `discriminator` [S]. The framework packages do not fit: flask-smorest (0.47.0) would mean rewriting `schemas.py` on marshmallow; apiflask (3.1.2) replaces the `Flask` class; flask-openapi3 (4.3.2) brings its own app classes; spectree (3.0.0) touches every route; flasgger (0.9.7.1) is unmaintained since 2023; flask-pydantic (0.14.0) writes no spec. schemathesis (4.30.0, Python >=3.10) is an optional test tool for the 3.12 CI job, since `test_fuzz_api_auth.py` already fuzzes with hypothesis. oasdiff is a Go binary, not on PyPI: `oasdiff breaking` lists changes that break clients and `changelog` lists all [S]; licence and release channel are [R].

Recommendation: no framework. A generator of about 150 lines (`planetgen/api/spec.py`) builds an OpenAPI 3.1 document for the remote contract only (roughly 12 to 15 routes: version, download, run create, id reservation, batch put, run status, seal, complete, recipes) from the Pydantic models. It adds no dependency that Python 3.9 cannot install. The 3.9 CI job compares the committed file by plain `json` equality; the 3.12 job runs `oasdiff`.

### 4.4 The rule (answer to API.4's open question)

1. `planetgen/api/version.py` holds `API_VERSION = "1.4"`, `MIN_CLIENT = "1.2"` and a history list (version, summary, converters).
2. `python -m planetgen.api.spec --check` regenerates the remote-contract spec and fails if it differs from the committed `docs/api/openapi.json`.
3. CI runs `oasdiff breaking` and `changelog` between PR base and head. A breaking diff requires a higher MAJOR; a non-breaking non-empty diff a higher MINOR; no diff an unchanged version. Description-only changes are excluded [S].
4. The check sits beside `scripts/bump_version.py --check-pr`. `stamp-version.yml` writes the release number into the compatibility table row after merge (unknown until then, so a hand-written column would always be wrong). No `changes/<name>.api.md` kind is needed: the spec diff is the record.
5. The docs compatibility table (API.4) is generated from the history list into `docs/api.md`. Columns: API version, release introduced, changed routes and payloads, client versions that work.

Without oasdiff, use the weaker rule "spec changed means version changed", with a reviewer checking major versus minor.

### 4.5 Version key and reproducibility (warning for API.17)

[reproducible-galaxies.md](reproducible-galaxies.md) section 4 puts the OS, architecture and Python micro version into the version key and says a seed reproduces a galaxy only on the exact setup. A remote run on a laptop with a different OS, architecture or Python never matches the server's key, even on an identical release, so API.17's check ("produces exactly what the server would make") holds only when the keys match. Design: `/api/version` returns the server's `version_key`; the client compares and says which parts differ; the server's sample re-run (API.8) is a hard check when the keys are equal and advisory otherwise. Whether a mismatch refuses the run is Boss's call. The sentence "No numpy or other numeric library is used" in section 4 of that note is out of date since the physics moved to scipy and astropy (GEN.66, `changes/physics-on-scipy-astropy.minor.md`), which makes the OS and architecture parts of the key matter more.

## 5. Key scopes (API.6, API.9)

### 5.1 Scopes

Prior art (GitHub fine-grained tokens, Stripe restricted keys `rk_...`, Google OAuth scopes [R]) points to `resource:action` plus a few coarse roles. This API has few resources, so the coarse set is enough:

| Scope | Allows | Implies |
|---|---|---|
| `read` | authenticated read routes and the key's own rate-limit bucket (reads are public anyway) | none |
| `generate` | submit recipes (API.18) and start generation jobs | `read` |
| `upload` | remote-run routes: download (API.12), reserve, batch upload, seal, status (API.3, API.14) | `read` |
| `admin` | every current write route (edits, jobs, admin stats, key listing) | all of the above |

A key's scopes are limited to what its owner's role may hold (Boss: "Only admin can upload"; user-level keys "access but not upload"). Today every owner is an admin, so the check is "creator is an admin" for `upload` and `generate`; USR.2 later maps roles to scopes. Account-management routes stay `session_only` whatever the scope ([../api.md](../api.md#authentication), TEST.44). `upload` is separate from `admin` so a key on a laptop cannot delete a sector or rename the galaxy if the laptop is lost. Refusal: `403 {"error":"this key lacks the 'upload' scope","required_scope":"upload"}`.

Open policy point: the galaxy seed (API.16, API.12) lets anyone with the code predict unfilled sectors. Default: those routes need `upload` or `generate`, not public.

### 5.2 Storage

- Keep SHA-256 of the full token: 256 bits of entropy make a slow hash pointless (bcrypt and scrypt are for low-entropy secrets) [R], and `UNIQUE(key_hash)` lookup leaks nothing because the attacker cannot choose the hash. If lookup is ever by id, use `hmac.compare_digest`.
- Keep `pg_` + 43 base64url characters; optionally append a 6-character CRC32 checksum (as GitHub's 2021 token format does [R]) so typos and junk are rejected without a database lookup. Old keys keep working.
- Add `key_prefix CHAR(8)` (shown in the key list and in log lines beside `key_id`) and `expires_at TIMESTAMP NULL`.
- Scopes in `admin_api_key_scopes (key_id BIGINT UNSIGNED, scope VARCHAR(32), PRIMARY KEY (key_id, scope), FOREIGN KEY ... ON DELETE CASCADE)`. Not a JSON column: MariaDB 11.4 stores JSON as LONGTEXT with a check while MySQL 8.4's is native, and control v10 (DB.13) just replaced a JSON text column with rows. A single `scope` column cannot say "upload but not admin" together with "generate". No CHECK on values, so a new scope needs no migration.
- Migration (control v11): create the table, add the columns, insert `admin` for every existing key.
- Rotation is create, switch, revoke. Expiry is chosen at creation (suggest 90 days for `upload` keys, none for `read`); the list shows "expires in N days"; an expired key is treated like a revoked one and logged `apikey.invalid`.
- Throttle `last_used_at`: `UPDATE ... WHERE key_hash=? AND (last_used_at IS NULL OR last_used_at < NOW() - INTERVAL 60 SECOND)`, and skip the commit when no row changed.

### 5.3 Per-key rate limits

Flask-Limiter runs its `key_func` in `before_request`, before `require_admin` runs, so the key is not yet known. Register a `before_request` hook ahead of `limiter.init_app(app)` that resolves the bearer token once (cached on `g`), and have `key_func` return `key:<id>` when present, else the IP. Never key by the raw header value (random strings would mint unlimited Redis keys). The per-IP default (200/day, 50/hour) must not apply to authenticated keys, or a `read` key is capped at 50 calls an hour. `cost` as a callable per limit [S] gives the bytes budget in 6.6. Consistent with [login-brute-force-protection.md](login-brute-force-protection.md): an API key never needs a TOTP code, key presentation is not under the login lockout, and a 256-bit key cannot be guessed, so the checksum pre-filter and the IP limit suffice. Flask-Limiter does not limit concurrency; a Redis counter with an expiry serves as the semaphore, and during an outage it falls back to per-process counting, acceptable with one process.

### 5.4 Checks for `tests/test_api_auth_sweep.py`

- `require_scope` sets `wrapped.scope_required`. A sweep asserts every `/api` route is exactly one of public (explicit allow-list), session-only or scoped, so a new route with none fails the build as a new unauthenticated write route fails TEST.43.
- A scope matrix from `app.url_map`: for each scope and route, 403 below the required scope and not 401/403 at or above it.
- Revoked, expired and owner-deleted keys give 401; an `upload` key cannot reach `session_only` routes; a `read` key never reaches a write route.
- No raw token in any column; key lists never include `key_hash`; a `read` or `upload` key gets its own limit bucket.

## 6. Upload limits and batches (API.3, API.7, API.12, API.13, API.14)

### 6.1 Measured sizes

Method [C]: `planetgen sector` on a throwaway MariaDB 10.11 (five sectors, 5 to 81 systems), every row exported with a foreign-key walk as JSON (one object per row, compact separators), then compressed. Also 40 systems generated in memory and measured as `StarSystem.to_dict()`.

| Sector | Systems | Rows | Raw JSON | gzip -6 | Raw per system |
|---|---|---|---|---|---|
| A | 5 | 310 | 387 KB | 84 KB | 77 KB |
| B | 8 | 351 | 416 KB | 84 KB | 52 KB |
| C | 20 | 1,167 | 1.38 MB | 285 KB | 69 KB |
| D | 72 | 4,170 | 4.82 MB | 979 KB | 67 KB |
| E | 81 | 4,409 | 5.15 MB | 1.04 MB | 64 KB |

In sector D, moons are 55% of the bytes (1,886 rows at 1.4 KB, about 26 per system), planets 22%, rogue planets 15%. Rows per system: moons 26, planets 9.8, rogue planets 9.2 (phenomena are a large share), stars 1.35, comets 0.8, belts 0.9. The server's own estimate for a 20-system sector is 1.1 MB in the database: 55 KB per system with indexes against 69 KB raw JSON.

In-memory `to_dict()` for 40 default-config systems: mean 58.6 KB raw, median 59.5 KB, p90 116 KB, max 155 KB; gzip -6 mean 13.8 KB (ratio 4.24); 0.04 s to generate and serialize one system.

Core sectors hold far more systems, so batch limits work by units, not by an assumed sector size. The row set leaves out population tables and caches the server rebuilds.

Recommendation: upload the `to_dict()` object graph in a small versioned envelope, not rows named after `schema.sql` columns. Column-shaped rows would force an API bump on nearly every schema migration (galaxy schema v65, control v10). The server saves the objects through the existing save code.

### 6.2 Compression (sector D, 4.82 MB, object-row JSON) [C]

| Codec | Size | Ratio | Compress MB/s | Decompress MB/s |
|---|---|---|---|---|
| gzip 1 | 1,161 KB | 4.15 | 31.5 | 65 |
| gzip 6 | 979 KB | 4.92 | 10.0 | 62 |
| gzip 9 | 955 KB | 5.04 | 7.5 | 98 |
| zstd 1 | 863 KB | 5.58 | 128 | 524 |
| zstd 19 | 735 KB | 6.56 | 0.2 | 1,311 |
| brotli 4 | 856 KB | 5.63 | 10.9 | 91 |
| brotli 11 | 703 KB | 6.85 | 0.1 | 87 |
| bz2 / xz | 740 / 743 KB | 6.5 | not timed | not timed |

One sample on a slow CPU. zstd 3 and 10 were also run (1,012 KB and 890 KB); zstd 3 coming out worse than zstd 1 is odd and should be re-run.

- Float-heavy data compresses 4 to 7 times, not the 10 to 20 times usual for JSON; the 17-digit doubles cannot be rounded without breaking exact reproduction.
- Columnar JSON halves the raw size (2.22 MB against 4.82 MB) but gzip ends only about 17% smaller (815 KB against 979 KB). Not worth a second format.
- gzip -6 is stdlib on every platform and Python 3.9 and understood by every proxy. zstd needs the `zstandard` package (0.25.0, BSD-3, Python >=3.9 [S]) or Python 3.14's stdlib module [R]; level 1 gains 12% over gzip -6 at 12 times the speed, which matters only if compression becomes the client's bottleneck. Ship gzip, define the envelope so `Content-Encoding: zstd` can be added without a version bump, skip bz2 and xz.
- Server work per 5 MB sector before validation: gunzip about 0.04 s, `json.loads` about 0.25 s (float-heavy), SHA-256 0.04 s [C]. Validation and inserts dominate and are measured when built.

### 6.3 Request decompression in Flask under mod_wsgi

- Flask and Werkzeug do not decode `Content-Encoding` on requests: a gzip body sent to `/api/auth/login` came back `400 request body must be a JSON object` [C].
- Apache can inflate with `SetInputFilter DEFLATE` (`DeflateInflateRatioLimit` default 200, burst 3 [S]); nginx and Caddy do not decompress requests natively [R]. Because deployments differ, decompress in the app so one path works everywhere.
- Implementation: accept only `Content-Encoding: gzip` (else `415`); require `Content-Length` (else `411`; chunked request bodies under mod_wsgi are unreliable [R]); read `request.stream` in 64 KiB pieces through `zlib.decompressobj(16 + zlib.MAX_WBITS)` with `decompress(chunk, max_length)`; answer `413` the moment the output passes the cap, or exceeds 40 times the input after the first MiB. A 300 MB zero payload (306 KB gzip) was cut at a 64 MiB cap in 0.94 s [C], hence the 32 MiB cap below. gzip cannot exceed about 1,032:1 (measured 1,028:1); zstd reached 32,676:1 and brotli 661,562:1 on zeros [C], so ratio guards matter more if those codecs are added. Tests: a bomb, a truncated stream, trailing garbage, a wrong encoding.
- Apply the cap per request: setting `request.max_content_length = 4 * 1024 * 1024` inside the upload view or a path-scoped `before_request` raised the limit for that route only (1 and 3 MB passed, 5 MB got 413, a 2 MiB route kept rejecting 3 MB) [C]. `_reject_oversized_body` compares to the app-wide config first, so it must become endpoint-aware or the override must run before it.
- Verify a `Content-Digest: sha-256=:...:` request header (RFC 9530 [R]) over the compressed bytes; it catches truncation and proxy damage before parsing.

### 6.4 Web server and proxy limits

| Layer | Setting | Default | Shipped | Needed for uploads |
|---|---|---|---|---|
| Apache | `LimitRequestBody` | 1 GiB since 2.4.54, unlimited before [S] | none | 5 MiB in `<Location "/api/uploads">` |
| nginx | `client_max_body_size` | 1m [R] | 1m | `location /api/uploads { client_max_body_size 5m; }` |
| Caddy | `request_body { max_size }` | none | 1MB | a path matcher with 5MB |
| Flask | `MAX_CONTENT_LENGTH` | none | 2 MiB | per-request override, 4 MiB |
| mod_wsgi / gunicorn | request timeout | none | 60 s | keep 60 s; size batches well inside it |

Apache's `LimitRequestBody` counts the compressed wire body unless the DEFLATE input filter is on. The 60 s timeout limits batches on slow links: 4 MiB at 1 Mbit/s takes 34 s and 8 MiB would not fit, so the client target is 1 to 2 MiB and it shrinks on a timeout. The five threads are shared with Server-Sent Event log streams, each holding a thread up to 40 seconds (section 10), so the concurrency limits below leave threads for pages.

### 6.5 Resumable batch design

tus 1.0.0 (HEAD for `Upload-Offset`, PATCH to append [S]), Google resumable uploads (256 KiB chunks [R]) and S3 multipart (numbered parts, list, complete or abort [R]) all resume one large object. Boss's units are many small, independent batches (15 to 1,000 KB), so the closest fit is the S3 model, with batch numbers and SHA-256 playing the part of ETags:

```
POST   /api/uploads                          create run (claims sectors, reserves ids)   Idempotency-Key
PUT    /api/uploads/<run>/batches/<seq>      one gzip batch; body hash in Content-Digest
GET    /api/uploads/<run>                    what the server holds: batch seqs + sha256, units sealed/finalized, id ranges left
POST   /api/uploads/<run>/units/<unit>/seal  manifest: row counts per table + unit SHA-256
POST   /api/uploads/<run>/complete           all units sealed; server finalizes the rest, reports
DELETE /api/uploads/<run>                    abort (admin clearing, ADM.13, uses the same call)
```

Re-sending a batch with the same seq and hash is a no-op returning the stored answer; the same seq with different bytes is `409`. The client keeps every batch on disk until `GET` lists its hash (API.13), so a crash, a lost response or a restarted server all resume by asking "what do you have". Verification and finalize of a sealed unit follow the queue rule (recommendation): `200` if it finishes inside the short wait, otherwise `202` with a job id on `GET /api/jobs/<id>`. Remote clients poll `GET /api/uploads/<run>`; the SSE log stream is a session-cookie page route.

Idempotency follows the IETF Idempotency-Key draft [S]: a missing key on a route that requires it is `400`; the same key with a different payload is `422`; a retry while the first is running is `409`. Use it on the non-batch POSTs; batch PUTs are idempotent by (run, seq, hash). Stored answers are kept 24 hours in Redis.

### 6.6 Numeric limit table (proposal for API.7)

| Limit | Value | Basis |
|---|---|---|
| Client batch target | 1 to 2 MiB compressed (about 5 to 10 MiB raw, 80 to 150 systems) | 60 s timeout, 1 Mbit/s uplink; 13.8 KB gzip per system [C] |
| Server hard limit, compressed body | 4 MiB, `413` above | 2 times the target; proxies set to 5 MiB on the upload path only |
| Server hard limit, decompressed | 32 MiB, and ratio over 40:1 after the first MiB | observed ratios 4 to 7 [C]; 306 KB bomb cut at cap [C] |
| JSON nesting depth | 12 on uploads, pre-scanned before `json.loads` | fixes the 500 on deep nesting [C] |
| Units per batch | whole units only; a larger unit spans batches and is sealed | a sector of 81 systems is about 1 MB gzip |
| Rows per request | 50,000 | sector E has 4,409 |
| Concurrent upload requests | 1 per key, 2 on the whole server | 5 threads in the shipped config; keep 3 for pages |
| Batch requests | 60 per minute per key | `WRITE_RATE_LIMIT` (10 per minute) would allow only 20 MiB per minute |
| Bytes budget | `cost = ceil(compressed_bytes / 256 KiB)`, limit 256 per minute (64 MiB per minute) per key | Flask-Limiter `cost` callable |
| Control calls (create, reserve, seal, status, complete) | 600 per hour per key | |
| Recipe submissions | 20 per minute per key, at most 5 queued jobs per key | 9.4 |
| Open runs | 2 per key | |
| Staging cap | 2 GiB on the spool volume and at most 10% of its free space; 1 GiB per run | only unsealed or unfinalized data; finalized units are deleted |
| Id reservation | per table, sized from the plan, 25% headroom, extendable by a reserve call | 7.2 |
| Stale flag (ADM.13) | no contact for 3 days: flag only, never auto-release | Boss: reserved ids stay until finished or an admin clears them |
| Idempotency answers kept | 24 h | |

Download side (API.12): `system_name_registry` rows are about 157 bytes of JSON, about 157 MB raw per million systems [C, extrapolated]. Serve registries paginated by cursor (`?after=<id>&limit=50000`), each page gzip'd, never one 60 s response.

### 6.7 What a client does at a limit

- `429` or `503`: wait `max(Retry-After, full-jitter backoff)`, a random value between 0 and min(60 s, 1 s times 2^attempt) (the AWS Architecture Blog scheme [R]). Trust `Retry-After` only on 429 and 503: Flask-Limiter sends `Retry-After: 3599` on a normal 200 [C]. Suppress it on non-429 responses (map `RATELIMIT_HEADER_RETRY_AFTER` or filter in `after_request`).
- `413` or a timeout: halve the batch (down to one unit) and resend only if nothing was committed; after three successes grow by 25% (AIMD).
- Other `5xx`: retry up to 6 times with the same batch (idempotent). Other `4xx`: stop and show the problem body. `409` on a claimed sector: skip it and report.
- Headers: Flask-Limiter emits `X-RateLimit-Limit`, `-Remaining`, `-Reset` (an epoch timestamp) and `Retry-After` [S][C]. The IETF draft replaces the trio with `RateLimit-Policy` and `RateLimit` structured fields [S] but is not an RFC and Flask-Limiter lacks it; keep the `X-` headers and treat all as hints.

## 7. Verification and staging (API.8, API.10, API.11, ADM.13)

### 7.1 Validation order (cheapest first, fail early)

1. Transport: authentication, scope, `Content-Length` against the limit, `Content-Encoding`, `Content-Digest`, concurrency slot.
2. Syntax: bounded decompress, depth pre-scan, `json.loads`. Failure is `400`.
3. Envelope (Pydantic): run id, seq, unit ids, `format_version`. An unknown version is `400` with the supported range.
4. Sequence: is (run, seq) already held? Same hash returns the stored answer; a different hash is `409`.
5. Reservation: the unit's sector is claimed by this run and row ids are inside its reserved ranges.
6. Per-unit schema (Pydantic models of the generator objects, `from_dict`): `422` with every error.
7. At seal time: counts per table and the unit hash match the manifest.
8. Server checks: names unique (needs the registry locks) and ADM.5's physics validator.
9. Finalize in one transaction through the normal save path; on any failure roll back and return the reasons.

### 7.2 Id blocks (API.10)

The existing allocator is already hi/lo: `_allocate_id` takes a block per process and table (64 doubling to 4,096) from `id_blocks` by `INSERT ... ON DUPLICATE KEY UPDATE next_id = LAST_INSERT_ID(GREATEST(next_id, floor) + size)` on its own autocommit connection, never below `MAX(id)+1`; rollbacks leave gaps. A remote run is the same pattern with a bigger lease: it asks for N ids per table (planned sectors times the rows per system in 6.1, plus 25%), gets contiguous ranges recorded in `upload_id_ranges (run_id, table_name, first_id, last_id)`, and asks for more with a `reserve` call. Unused ids are abandoned; gaps are normal. Alternatives rejected: ULID, UUIDv7 or Snowflake ids (BIGINT columns and every table's ordering would be rethought) and client-local ids remapped at finalize (Boss decided the server reserves blocks, and population streams still key on row ids, GEN.57).

Name clashes: registries downloaded once go stale while other uploads land. The server renames the system on finalize and its planets and moons through the admin rename code, and returns a rename map so the client's cache stays consistent.

### 7.3 Staging (API.11)

A spool directory (config `uploads.dir`, like `jobs.dir`) holds received gzip batches as files named by run and seq. Control-DB tables: `upload_runs` (id, owner, key_id, galaxy db, state, created, last_contact), `upload_batches` (run, seq, sha256, bytes, received_at), `upload_units` (run, unit key, state open, sealed, verified, finalized or rejected, manifest hash, errors), `upload_claims` (run, sector address, state) and `upload_id_ranges`. A unit's JSON is assembled in memory (a few MB) at seal. Mirror staging tables in the galaxy schema, as API.11's text says, are not recommended: about 30 content tables and 65 migrations would each need a staging twin with the foreign-key and trigger behaviour again. Boss's requirement is "write to the database in a controlled way only complete units", so the unit boundary, not row-level staging, is what matters. A unit's batches are deleted after it finalizes and its manifest kept. These are control v12, separate from the scope migration.

### 7.4 Hash manifests and API.8's open question

Per batch: SHA-256 of the received bytes. Per unit: SHA-256 over the canonical JSON plus row counts per table. Per run: the sorted list of unit hashes hashed at `complete`. A Merkle tree pays off only to prove one leaf among millions; with tens of thousands of units a flat list is simpler. For API.17 the server may re-run a random sample of sectors and compare unit hashes, valid only when version keys match (4.5).

API.8's open question, default proposal: the server corrects only what is deterministic and cosmetic, a clashing name (rename plus dependents) and derived columns a pure function recomputes. It rejects the unit for anything else (outside the reserved sector or id range, missing children, failed physics checks, wrong format version) with reasons, and the client regenerates it. A rejection does not release the claim.

## 8. Logging every API call (API.15)

### 8.1 Fields and hook

Time (UTC, ms), request id (also `X-Request-Id` and in error bodies), method, matched route (`request.url_rule.rule`, `-` for a 404) plus the raw path truncated to 200 characters, status, duration ms, bytes in and out, user (key owner, signed-in admin, `god` or `-`), `via` (`key`, `session`, `anonymous`, `console`), `key_id` and `key_prefix` (never the token), `internal=1` for page-to-API calls made in-process, client IP, user agent (120 characters), the client's API version header, and `ratelimited=1`. Boss's required parts (time, route, user, how it came in, status) are all there. The docs say status codes are HTTP codes for parity with the web server's logs, and that requests the web server rejects before Flask (a 413 from `LimitRequestBody`, TLS failures) appear only in its log.

"Console": `planetgen` commands write to the database directly and never call the HTTP API, so nothing on a request path is `god` today. Default: `god` lines come from API-shaped actions run by CLI tools or an in-process client with no request context; the existing `GEN` lines stay.

Hook `after_request`: it sees every response including error-handler ones (a probe recorded 404 with endpoint `None`, 401, 413, 400 and 429, the 429 carrying `X-RateLimit-*`) [C]. `teardown_request` is only for crashes with no response. WSGI middleware loses `g` (user, key).

In-process calls: opening `/`, `/sectors`, `/sector/1`, `/system/1`, `/search` and `/galaxy` each made 6 in-process API calls (3 `galaxy/changes`, `population`, `species`, the page's data call) [C; the pages returned 404 on that database, so real pages may make more]. Logged as ordinary lines they multiply volume about six times with lines that look like traffic. Mark them `internal=1` and add a setting `api_log.internal` (`all` by default, per Boss's "standard will have all API calls logged"; `off` drops them).

### 8.2 Cost and storage

Measured [C]: `activity_log.event` takes 78 to 126 microseconds per call (8,000 to 13,000 per second); one autocommit MariaDB INSERT 5.9 ms (slow sandbox disk; expect 0.5 to 2 ms on an SSD); 500 rows per commit 0.43 ms per row; about 98 bytes per row with two secondary indexes. At 1 million calls a day that is about 100 MB a day and 3.6 GB a year.

1. **Now (API.15):** add `API` to `CATEGORIES` (and its row in `docs/config.md`) and write one line per call to the activity log, with its rotation and fail2ban-compatible format and the rule that writing never stops the program. Add `X-Request-Id`. No schema change, no user accounts needed.
2. **When an admin page needs search:** an `api_calls` table in the control DB (same columns) filled by a background thread draining a bounded `queue.Queue`, 500 rows per commit; when full, drop and count. The deployment is one process with threads, so an in-process queue is sound. Redis would add a failure mode for no gain; synchronous inserts add milliseconds to every call and stack up across 5 threads.

Retention: full IPs 30 days, then truncated (IPv4 /24, IPv6 /48), rows kept 13 months; purge by batched `DELETE ... LIMIT 10000` from the daily maintenance run (OPS.16). Partition by month only past tens of millions of rows. Privacy: an IP address is personal data under GDPR where the operator can link it to a person, and dynamic addresses count (CJEU, Breyer, C-582/14 [R]); no statutory retention number exists for security logs, and 30 days full plus truncation is a defensible default for a hobby site. The activity log already holds IPs under logrotate (100 MB, 30 copies). State this in the docs. Never log request bodies, tokens or `Authorization` headers.

## 9. Recipes (API.18, API.19)

### 9.1 Pydantic and the shape

Pydantic v2 models are the source of truth, with JSON Schema exported by `model_json_schema` for the docs and the OpenAPI file. The repo already depends on Pydantic >=2.7, validates every body with it (`api/schemas.py`) and the Generate page uses `TypeAdapter`. Python 3.9 works with 2.13.5 [S] under two rules: write unions as `Union[...]` and `Optional[...]` (`X | Y` fails at runtime on 3.9), and stay within 2.13.x features. The `--system-file` format ([../system-file-format.md](../system-file-format.md): tri-state options, per-orbit `slots`) is already a recipe; the schema starts from it and adds `fixed / range / random` for any numeric field. Prototype [C]:

```python
Dim = Annotated[Union[Fixed, Range, Random], Field(discriminator="mode")]   # mode: fixed|range|random
class Planet(BaseModel): kind: Literal["planet"]; mass_earth: Optional[Dim]; moons: List[Moon] = Field(max_length=100)
class System(BaseModel): kind: Literal["system"]; recipe_version: int; star_type: Optional[str];
                          slots: List[Annotated[Union[Planet, Belt], Field(discriminator="kind")]] = Field(max_length=500)
```

Top level: `{"recipe_version": 1, "kind": "sector|system|planet|moon|phenomenon|galaxy", ...}`, with `phenomenon` carrying a `type` discriminator. `Range` has `distribution` (`uniform`, `log-uniform`, `normal`) and a check `min <= max`. Models use `extra="forbid"` and `strict=True`, as `schemas.py` does. Pydantic reports a bad tag as `union_tag_invalid` with the allowed list and reports all errors at once with `loc` paths such as `slots.0.planet.mass_earth.range` [C]. `recipe_version` is separate from the API version; older recipe versions go through `from_v1(...)` upgraders. Galaxy-scale recipes (API.19) are lists of region recipes applied one at a time through the same endpoint, not one giant document.

### 9.2 Determinism and the seed scheme

```
unit seed = SHA-256(galaxy seed || "recipe:" || recipe_hash || ":" || seed)
```

`recipe_hash` is the SHA-256 of the canonical JSON of the validated recipe with defaults filled in (`model_dump(mode="json", by_alias=True)`, sorted keys, compact separators, no NaN), so omitted and explicit defaults hash the same. `seed` is the recipe's optional explicit seed or a drawn one returned in the result (`seed_used`), as `planetgen plan` draws a galaxy seed. Resubmitting with `seed_used` on the same release reproduces the output; the log line states recipe hash, seed and version key. Float `repr` in `json.dumps` is the same shortest round-trip form on 3.9 to 3.13 [R]; pin it with a golden-hash test. This is a unit seed as in [reproducible-galaxies.md](reproducible-galaxies.md) section 3: the recipe job draws through `planetgen/util/draw.py` under `seeded`. Limitation: the unit stream is sequential, so any recipe edit re-rolls everything not fixed (per-field substreams would need GEN.56's helpers keyed by path); to change one field and keep the rest, fix the others.

### 9.3 400 versus 422

RFC 9110 section 15.5.21 defines `422 Unprocessable Content`: content type and syntax are understood but the instructions cannot be processed [S via the Idempotency draft; wording R]. Two layers: an envelope check gives 400, then the kind's model plus physics checks give 422.

- **400, gibberish:** not valid UTF-8 or JSON; not an object; nested deeper than the cap (including the `RecursionError` case); no `recipe_version` or `kind`, or an unknown `kind`.
- **422, validation:** a known kind with a wrong field: unknown field, wrong type, out of range, `min > max`, a list too long, a physics or consistency rule broken (`+habitable_world` with an O star), or a generation that cannot satisfy the recipe after its retries (the CLI "tried five times" case).
- **413** size cap, **415** wrong `Content-Type`, **429** rate limits.

Body: RFC 9457 `application/problem+json` (obsoletes RFC 7807 [R]) with `type`, `title`, `status`, `detail`, and extension members `errors` (each with `pointer`, a JSON Pointer such as `/slots/3/mass_earth`, `code` and `message`), `log` (the last 100 lines or 16 KB of the command's log), `request_id`, and `truncated` when more than 50 errors exist. The existing `error` and `errors:[{field,message}]` members stay so current clients keep working: `ApiError` gets optional `type` and `log`, and a `ValidationFailed` subclass supplies 422. The log is filtered to generation log lines (no file paths, SQL or tracebacks), consistent with the 500 handler never leaking internals. Existing routes keep `400` for body validation (changing them is an API-major event under 4.4), and [../api.md](../api.md#errors) should state the split.

### 9.4 Cost limits

Validation cost [C]: 0.1 ms for a 3 KB system recipe, 24.5 ms for 91 KB, 2.35 s for a 4 MB recipe with 500 slots times 100 moons (about 0.6 ms per KB). Proposed caps: body 256 KiB (about 70 ms), depth 12, at most 500 slots (`tuning.ABSOLUTE_MAX_SYSTEM_OBJECTS`), 100 moons per planet, sector and neighbourhood recipes bounded by `generation/limits.py` (`MAX_GENERATE_RADIUS_PC = 200`, `MAX_GENERATE_LIMIT`), 5 queued recipe jobs per key, 50 reported errors. Recipes run as queue jobs, so the request thread only validates and enqueues; a slow one answers `202` with `GET /api/jobs/<id>`.

### 9.5 Bug found while testing

`require_json_body` uses `request.get_json(silent=True)`, which swallows `ValueError` but not `RecursionError`. A 100,000-deep `[` body (200 KB, under the 2 MiB limit) to `/api/auth/login` returned `500 {"error":"internal server error"}` with a traceback; depths 10, 500 and 5,000 gave 400 [C]. Fix: pre-scan nesting depth in a byte loop (or catch `RecursionError`) and answer 400. It matters more for recipes, which expect nesting.

## 10. The Web UX and Job Management Guide: what carries over

The guide ("Web UX and Job Management Guide.md", titled "Single-Server Migration Guide & Codebase Modernization Outline") was checked against the code and [library-migration.md](library-migration.md).

- **Superseded:** the broker-less `ProcessPoolExecutor` with SQLite state, `future.cancel()` and the SQLite WAL mitigation (Phase 1, "How should single-server queue management be structured"). Boss chose Redis; jobs are RQ jobs (PERF.24, built for the routes in [../api.md](../api.md#queued-jobs)) and cancel and pause use RQ's stop command.
- **Built with changes:** SSE log streaming (Phase 2, ADM.22, ADM.40). The stream is a page route, `/admin/generate/jobs/<id>/stream`, not `/api/jobs/<id>/stream`; events carry the byte offset as their id so a reconnect resumes; the log is a file read by offset, not an in-memory queue per job id, so it works across Apache's and the workers' processes; the stream closes itself every 40 s to stay under the daemon timeout. API clients get `GET /api/jobs/<id>` (state, result, error), not a stream.
- **Holds:** WSGI worker starvation. Under Apache one open stream holds one of the 5 threads; the deployment docs set threaded mod_wsgi daemon processes, and the upload concurrency limits above leave threads for pages and streams.
- **Not API:** Phases 3 and 4 (TanStack tables, Shoelace, 3D tiles, camera-relative coordinates) are UI work, built (UX.40, UX.41, MAP.102).

## 11. Evidence notes and sources

[C] results come from research scripts that are not in the repo: sector sizes and compression, per-system dict sizes, activity-log and INSERT timings, in-process call counts, deep JSON and size-limit behaviour, `after_request` coverage, `Retry-After` on 200, recipe validation timing, the per-request `max_content_length` override.

[R] claims to verify when paper and RFC access is allowed:

- RFC 9745 (Deprecation) number, date and format; RFC 8594 (Sunset); the `Link` relations; RFC 9530 `Content-Digest`; RFC 9457 obsoleting 7807; RFC 6648 on `X-`; RFC 9110 section 15.5.21 text for 422 and the 426 semantics.
- Stripe, GitHub (fine-grained tokens, checksum token format) and Google resumable-upload details; S3 multipart limits (5 MiB parts, 10,000 parts).
- Defaults: nginx `client_max_body_size 1m`; whether nginx and Caddy lack request decompression; chunked request bodies under mod_wsgi.
- AWS full-jitter backoff, CJEU Breyer C-582/14, any regulator guidance on log retention.
- oasdiff licence and release channel; schemathesis behaviour on Flask.
- Flask-Limiter 3.11.0 (the Python 3.9 pin) against 4.1.1 for `cost`, `key_func` and `Retry-After`; only 4.1.1 was tested.
- The IETF ratelimit draft's current number; JSON float representation identical across Python 3.9 to 3.13; a re-run of the zstd 3 sample.

Sources:

- tus protocol 1.0.0: https://raw.githubusercontent.com/tus/tus-resumable-upload-protocol/main/protocol.md
- Apache mod_deflate: https://raw.githubusercontent.com/apache/httpd/trunk/docs/manual/mod/mod_deflate.xml
- Apache core (LimitRequestBody default): https://raw.githubusercontent.com/apache/httpd/trunk/docs/manual/mod/core.xml
- IETF RateLimit header fields draft: https://raw.githubusercontent.com/ietf-wg-httpapi/ratelimit-headers/main/draft-ietf-httpapi-ratelimit-headers.md
- IETF Idempotency-Key draft: https://raw.githubusercontent.com/ietf-wg-httpapi/idempotency/main/draft-ietf-httpapi-idempotency-key-header.md
- Pydantic unions and discriminators: https://raw.githubusercontent.com/pydantic/pydantic/main/docs/concepts/unions.md
- oasdiff breaking changes: https://raw.githubusercontent.com/oasdiff/oasdiff/main/docs/BREAKING-CHANGES.md
- PyPI JSON (https://pypi.org/pypi/<name>/json, read 2026-10-09) for pydantic, flask, flask-limiter, limits, flask-smorest, apiflask, flasgger, spectree, flask-pydantic, flask-openapi3, schemathesis, zstandard, brotli; oasdiff returned 404 on PyPI.
- Flask-Limiter 4.1.1 source as installed (`constants.py`, `_extension.py`).
