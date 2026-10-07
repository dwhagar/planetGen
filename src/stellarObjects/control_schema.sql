-- stellarObjects/control_schema.sql
--
-- MySQL (InnoDB) schema for planetGen's *control plane* -- admin
-- identities, web sessions, API keys, and the audit trail of every
-- write/admin action -- kept in one dedicated schema
-- (`PLANETGEN_CONTROL_DATABASE`, default `planetgen_control`) separate
-- from every per-galaxy content schema `schema.sql` describes.
--
-- Why separate: a deployment can host several galaxy databases sharing
-- one MySQL server (`stellarObjects._db.list_databases`, prefix-filtered
-- -- one schema per campaign/galaxy). Admin accounts describe *who can
-- administer this deployment*, not *who owns one galaxy*, so duplicating
-- `admin_users` into every content schema would mean re-registering the
-- same handful of admins in each one and having no single source of
-- truth for "is this session/API key still valid" across them. This
-- schema is applied once, independent of how many content schemas exist.
--
-- Versioned the same way `schema.sql` is (see that file's header and
-- `stellarObjects/_db.py`'s `SCHEMA_VERSION`/`migrate_database`), but with
-- its own independent version counter (`control_schema_migrations`) --
-- this schema's shape has nothing to do with a galaxy's own generated
-- content, so there is no reason its version number should ever need to
-- track `schema.sql`'s.
--
-- Secrets are never stored in plaintext or in a reversible form:
--   - `admin_users.password_hash` -- `werkzeug.security.
--     generate_password_hash` (PBKDF2/scrypt, salted), verified via
--     `check_password_hash`. See `planetgen/admin/auth.py`.
--   - `admin_sessions.token_hash` / `admin_api_keys.key_hash` -- SHA-256
--     hex digest of a `secrets.token_urlsafe` value. The raw token/key is
--     shown to the caller exactly once (at login / at creation) and never
--     persisted -- only its hash is stored, so a database read alone
--     can't be used to impersonate a session or API key, the same
--     "store a hash, not the secret" principle applied to the password
--     above.
CREATE TABLE IF NOT EXISTS control_schema_migrations (
    version      INT NOT NULL PRIMARY KEY,
    applied_at   TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- One row per admin. Deliberately no roles/permissions column -- per the
-- brief, this project only ever has a handful of admins, all with the
-- same full access, not a general user-accounts system.
--
-- `must_change_credentials` starts TRUE on the seeded `admin`/`password`
-- bootstrap row (see `adminAuth.ensure_default_admin`) and on any admin
-- created afterward with a caller-chosen initial password -- cleared only
-- by a successful `POST /api/auth/change-credentials`. While TRUE for the
-- calling admin, every write/admin endpoint except change-credentials
-- itself refuses the request (see `html/api/auth.py`'s
-- `require_fresh_credentials`), not just a UI reminder.
CREATE TABLE IF NOT EXISTS admin_users (
    id                        BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    username                  VARCHAR(64) NOT NULL,
    password_hash             VARCHAR(255) NOT NULL,
    must_change_credentials   TINYINT(1) NOT NULL DEFAULT 1 CHECK (must_change_credentials IN (0, 1)),
    created_at                TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    last_login_at             TIMESTAMP NULL,

    UNIQUE (username)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Web-UI login sessions (the `HttpOnly`/`Secure`/`SameSite=Strict` cookie
-- `POST /api/auth/login` sets). `expires_at` is a fixed lifetime set at
-- creation (no sliding-window renewal), so a stolen cookie has a bounded
-- useful life regardless of ongoing use -- deliberately simple over a
-- refresh-token scheme for a "a few admins" deployment.
CREATE TABLE IF NOT EXISTS admin_sessions (
    id               BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    admin_user_id    BIGINT UNSIGNED NOT NULL,
    token_hash       CHAR(64) NOT NULL,
    created_at       TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    expires_at       TIMESTAMP NOT NULL,
    last_seen_at     TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT fk_admin_sessions_admin_user
        FOREIGN KEY (admin_user_id) REFERENCES admin_users(id) ON DELETE CASCADE,
    UNIQUE (token_hash),
    KEY idx_admin_sessions_admin_user_id (admin_user_id),
    KEY idx_admin_sessions_expires_at (expires_at)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- API keys for programmatic (non-browser) callers -- `Authorization:
-- Bearer <key>`. `revoked_at` (not a DELETE) keeps a revoked key's audit
-- trail (`admin_audit_log.detail` can reference it, and `last_used_at`
-- stays visible in `GET /api/auth/api-keys` after revocation) instead of
-- losing history the moment a key is retired.
CREATE TABLE IF NOT EXISTS admin_api_keys (
    id               BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    admin_user_id    BIGINT UNSIGNED NOT NULL,
    label            VARCHAR(128) NOT NULL,
    key_hash         CHAR(64) NOT NULL,
    created_at       TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    last_used_at     TIMESTAMP NULL,
    revoked_at       TIMESTAMP NULL,

    CONSTRAINT fk_admin_api_keys_admin_user
        FOREIGN KEY (admin_user_id) REFERENCES admin_users(id) ON DELETE CASCADE,
    UNIQUE (key_hash),
    KEY idx_admin_api_keys_admin_user_id (admin_user_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Every write/admin action, one row each -- written by the same request
-- that performs the action (`planetgen.admin.auth.record_audit`), not
-- derived after the fact from `schema.sql`'s own tables (which don't
-- carry a "who changed this" column). `admin_username` is captured at
-- write time (denormalized, alongside `admin_user_id`) so the trail stays
-- readable even if that admin is later renamed or removed --
-- `admin_user_id` alone would go stale (`ON DELETE SET NULL`) in exactly
-- the case an audit log matters most.
CREATE TABLE IF NOT EXISTS admin_audit_log (
    id               BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    admin_user_id    BIGINT UNSIGNED,
    admin_username   VARCHAR(64) NOT NULL,
    action            VARCHAR(64) NOT NULL,   -- e.g. "sector.create", "system.delete"
    target            VARCHAR(255),           -- e.g. "sector:42", "system:1001"
    detail            TEXT,
    created_at        TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT fk_admin_audit_log_admin_user
        FOREIGN KEY (admin_user_id) REFERENCES admin_users(id) ON DELETE SET NULL,
    KEY idx_admin_audit_log_admin_user_id (admin_user_id),
    KEY idx_admin_audit_log_created_at (created_at)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v2 (SEC.1, SEC.21): failed-login counters and lockouts, one row per
-- client address (`scope` 'ip'; an IPv6 address by its /64) or username
-- (`scope` 'user', case-folded), shared by every worker process and kept
-- across restarts. See `planetgen/admin/throttle.py` for the rules.
-- Times are Unix seconds (DOUBLE), compared with the web server's own
-- clock; 0 means never. Idle rows are deleted after a week.
CREATE TABLE IF NOT EXISTS login_throttle (
    scope             VARCHAR(8) NOT NULL,
    subject           VARCHAR(128) NOT NULL,
    failures          INT UNSIGNED NOT NULL DEFAULT 0,
    level             INT UNSIGNED NOT NULL DEFAULT 0,  -- lockouts so far (doubles the next one)
    locked_until      DOUBLE NOT NULL DEFAULT 0,
    last_failure_at   DOUBLE NOT NULL DEFAULT 0,
    last_lockout_at   DOUBLE NOT NULL DEFAULT 0,

    PRIMARY KEY (scope, subject),
    KEY idx_login_throttle_locked_until (locked_until),
    KEY idx_login_throttle_last_failure_at (last_failure_at)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v3 (SEC.22): trusted devices. A successful login gives the browser a
-- long-lived device cookie (only its SHA-256 hash stored here); a login
-- from a browser holding a valid one for that username skips the
-- per-username lock, so someone failing on purpose can't lock the real
-- admin out (the per-address lock still applies). Changing credentials
-- deletes the admin's rows.
CREATE TABLE IF NOT EXISTS admin_devices (
    id               BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    admin_user_id    BIGINT UNSIGNED NOT NULL,
    token_hash       CHAR(64) NOT NULL,
    created_at       TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    expires_at       TIMESTAMP NOT NULL,
    last_used_at     TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT fk_admin_devices_admin_user
        FOREIGN KEY (admin_user_id) REFERENCES admin_users(id) ON DELETE CASCADE,
    UNIQUE (token_hash),
    KEY idx_admin_devices_admin_user_id (admin_user_id),
    KEY idx_admin_devices_expires_at (expires_at)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v4 (SEC.26): optional two-factor sign-in. One row per admin who has
-- started setting up an authenticator app; `enabled_at` is NULL until the
-- first code confirms it. `secret` is the base32 TOTP key itself (an
-- authenticator code can't be checked against a hash of it), so this
-- table is as sensitive as the database password. `last_step` is the
-- newest 30-second step a code was accepted for; older or equal ones are
-- refused, so no code works twice.
CREATE TABLE IF NOT EXISTS admin_totp (
    admin_user_id    BIGINT UNSIGNED NOT NULL PRIMARY KEY,
    secret           VARCHAR(64) NOT NULL,
    created_at       TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    enabled_at       TIMESTAMP NULL,
    last_step        BIGINT NOT NULL DEFAULT 0,

    CONSTRAINT fk_admin_totp_admin_user
        FOREIGN KEY (admin_user_id) REFERENCES admin_users(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v4 (SEC.26): single-use recovery codes, shown once when two-factor
-- sign-in is turned on; only their SHA-256 hashes are kept.
CREATE TABLE IF NOT EXISTS admin_recovery_codes (
    id               BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    admin_user_id    BIGINT UNSIGNED NOT NULL,
    code_hash        CHAR(64) NOT NULL,
    created_at       TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    used_at          TIMESTAMP NULL,

    CONSTRAINT fk_admin_recovery_codes_admin_user
        FOREIGN KEY (admin_user_id) REFERENCES admin_users(id) ON DELETE CASCADE,
    UNIQUE (code_hash),
    KEY idx_admin_recovery_codes_admin_user_id (admin_user_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v5 (PERF.8): the generation work queue (`stellarObjects/workQueue.py`).
-- `work_jobs` is one row per run that queued work (a `generate.py`
-- galaxy, sector or plan run, from the command line or the Generate
-- page); `work_tasks` is one row per unit it handed to the worker pool
-- (a sector to fill, a layer of bright stars to scatter), with its
-- timings and result; `work_lease` is the single row saying which run's
-- pool is using the machine's CPUs right now. The run holding the lease
-- refreshes `heartbeat_at` every few seconds and clears `holder` when it
-- finishes; a lease or job whose heartbeat is older than 30 s belongs to
-- a run that died, and the next run takes the lease over and marks that
-- run's unfinished tasks `cancelled`. Old jobs and their tasks are
-- deleted after a week.
CREATE TABLE IF NOT EXISTS work_jobs (
    id               VARCHAR(32) NOT NULL PRIMARY KEY,
    title            VARCHAR(255) NOT NULL,
    holder           VARCHAR(255) NOT NULL,   -- host:pid:token of the run
    state            VARCHAR(16) NOT NULL,    -- waiting, running, paused, done, failed, cancelled
    workers          INT UNSIGNED NOT NULL,   -- pool size, 0 for a node with no pool of its own
    tasks_queued     INT UNSIGNED NOT NULL DEFAULT 0,
    tasks_done       INT UNSIGNED NOT NULL DEFAULT 0,
    tasks_failed     INT UNSIGNED NOT NULL DEFAULT 0,
    created_at       DATETIME(6) NOT NULL,
    started_at       DATETIME(6) NULL,
    finished_at      DATETIME(6) NULL,
    heartbeat_at     DATETIME(6) NOT NULL,
    -- v7 (ADM.12): the job tree, see the comment after this table.
    parent_id        VARCHAR(32) NULL,
    root_id          VARCHAR(32) NULL,
    kind             VARCHAR(32) NOT NULL DEFAULT 'queue',
    seconds          DOUBLE NULL,             -- finished_at - started_at
    tasks_total      INT UNSIGNED NULL,       -- tasks the run said it would queue
    web_job_id       VARCHAR(32) NULL,        -- the Generate page's job (web/jobs.py)
    database_name    VARCHAR(64) NULL,        -- the galaxy database it wrote
    argv             TEXT NULL,               -- JSON command line, without --mysql-* options
    control          VARCHAR(16) NULL,        -- ADM.10: 'pause' or 'cancel', asked from the web

    KEY idx_work_jobs_state (state),
    KEY idx_work_jobs_created_at (created_at),
    KEY idx_work_jobs_root (root_id),
    CONSTRAINT fk_work_jobs_parent FOREIGN KEY (parent_id) REFERENCES work_jobs(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v7 (ADM.12): every job is a tree. A `generate.py` run is the root
-- (`kind` = its command), or, when the Generate page started it, a
-- child of that page job's "step" node, whose parent is the page job
-- itself (`kind` 'web-job', `web_job_id` its directory). Phases of a run
-- (the skeleton, the bright stars, population) are nodes under it, and
-- each work queue (`kind` 'queue') is a node whose leaves are its
-- `work_tasks` rows (a sector, a layer of bright stars). Every node keeps
-- its own start, end and `seconds`; a parent's totals are added up from
-- its children when the tree is read (`workQueue.load_tree`), so writing
-- a task touches only its own queue's row. `root_id` fetches a whole
-- tree at once; deleting a root deletes its subtree and tasks. `control`
-- is how the admin queue page (ADM.10) asks a node to pause or cancel.

CREATE TABLE IF NOT EXISTS work_tasks (
    id               BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    job_id           VARCHAR(32) NOT NULL,
    kind             VARCHAR(32) NOT NULL,    -- e.g. "sector"
    task_key         VARCHAR(255) NOT NULL,   -- e.g. "12,0,345" (ring, layer, slot)
    weight           DOUBLE NOT NULL DEFAULT 1,
    state            VARCHAR(16) NOT NULL,    -- queued, running, done, failed, cancelled
    created_at       DATETIME(6) NOT NULL,
    started_at       DATETIME(6) NULL,
    finished_at      DATETIME(6) NULL,
    seconds          DOUBLE NULL,             -- time the worker spent on it
    result           TEXT NULL,               -- short JSON summary
    error            TEXT NULL,

    CONSTRAINT fk_work_tasks_job FOREIGN KEY (job_id) REFERENCES work_jobs(id) ON DELETE CASCADE,
    KEY idx_work_tasks_job_state (job_id, state)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS work_lease (
    id               TINYINT UNSIGNED NOT NULL PRIMARY KEY,   -- always 1
    holder           VARCHAR(255) NULL,       -- NULL: free
    job_id           VARCHAR(32) NULL,
    heartbeat_at     DATETIME(6) NULL,
    -- v7 (ADM.10): "Pause the queue": while set, no run takes the lease
    -- or hands out a task, until an admin resumes the queue.
    paused           TINYINT(1) NOT NULL DEFAULT 0,
    paused_by        VARCHAR(64) NULL,
    paused_at        DATETIME(6) NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v6 (PERF.3, PERF.10): how fast this server generates and how much
-- space a galaxy takes (`planetgen/generation/stats.py`).
-- `generation_stats` is one row per kind of task ("sector" fill, or a
-- "scatter" layer of bright stars) and log-scale density bucket (two per
-- decade from 0.01, open-ended upward), each a decaying average over
-- every task that ever finished in it. `generation_size` is one row per
-- galaxy database: its bytes per star system, measured from its own
-- tables after each run. Bulk generation reads both for its size and
-- time estimate. A galaxy reset keeps them (they describe the server).
CREATE TABLE IF NOT EXISTS generation_stats (
    kind                 VARCHAR(16) NOT NULL,     -- sector, scatter
    bucket               INT NOT NULL,             -- floor(2 * log10(density / 0.01))
    density_low          DOUBLE NOT NULL,
    density_high         DOUBLE NOT NULL,
    samples              BIGINT UNSIGNED NOT NULL DEFAULT 0,
    seconds_per_task     DOUBLE NOT NULL DEFAULT 0,  -- worker wall time
    seconds_per_system   DOUBLE NOT NULL DEFAULT 0,
    systems_per_task     DOUBLE NOT NULL DEFAULT 0,
    stars_per_system     DOUBLE NOT NULL DEFAULT 0,
    max_density          DOUBLE NOT NULL DEFAULT 0,
    updated_at           DATETIME(6) NOT NULL,

    PRIMARY KEY (kind, bucket)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS generation_size (
    database_name        VARCHAR(64) NOT NULL PRIMARY KEY,
    bytes_per_system     DOUBLE NOT NULL,
    systems              BIGINT UNSIGNED NOT NULL,
    total_bytes          BIGINT UNSIGNED NOT NULL,
    measured_at          DATETIME(6) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
