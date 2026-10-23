-- planetgen/db/control_schema.sql
--
-- MySQL (InnoDB) schema for planetGen's *control plane* -- admin
-- identities, web sessions, API keys, and the audit trail of every
-- write/admin action -- kept in one dedicated schema
-- (`PLANETGEN_CONTROL_DATABASE`, default `planetgen_control`) separate
-- from every per-galaxy content schema `schema.sql` describes.
--
-- Why separate: a deployment can host several galaxy databases sharing
-- one MySQL server (`planetgen.db.store.list_databases`, prefix-filtered
-- -- one schema per campaign/galaxy). Admin accounts describe *who can
-- administer this deployment*, not *who owns one galaxy*, so duplicating
-- `admin_users` into every content schema would mean re-registering the
-- same handful of admins in each one and having no single source of
-- truth for "is this session/API key still valid" across them. This
-- schema is applied once, independent of how many content schemas exist.
--
-- Versioned the same way `schema.sql` is (see that file's header and
-- `planetgen/db/store.py`'s `SCHEMA_VERSION`/`migrate_database`), but with
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

-- One row per account. Until v15 every row was an admin with the same full
-- access; v15 (DB.58, for the accounts of USR.2 and the roles of 10.0) adds
-- `role` (`user`, `moderator`, `admin` or `owner`), `email` (NULL, unique) and
-- `status` (`active` or `disabled`). A moderator is accepted by the column's
-- CHECK and has no permissions until 10.0. There is exactly one `owner`: the
-- generated column `owner_flag` is 1 for an owner row and NULL for every
-- other, and its UNIQUE key lets only one row hold the 1. The row the install
-- seeds is the owner; an existing control database's lowest account id became
-- the owner when v15 was applied (`store._promote_control_owner`).
--
-- `must_change_credentials` starts TRUE on the seeded `admin`/`password`
-- bootstrap row (see `adminAuth.ensure_default_admin`) and on any admin
-- created afterward with a caller-chosen initial password -- cleared only
-- by a successful `POST /api/auth/change-credentials`. While TRUE for the
-- calling admin, every write/admin endpoint except change-credentials
-- itself refuses the request (see `planetgen/api/auth.py`'s
-- `require_fresh_credentials`), not just a UI reminder.
CREATE TABLE IF NOT EXISTS admin_users (
    id                        BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    username                  VARCHAR(64) NOT NULL,
    password_hash             VARCHAR(255) NOT NULL,
    must_change_credentials   TINYINT(1) NOT NULL DEFAULT 1 CHECK (must_change_credentials IN (0, 1)),
    created_at                TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    last_login_at             TIMESTAMP NULL,
    email                     VARCHAR(254) NULL,
    role                      VARCHAR(16) NOT NULL DEFAULT 'admin',
    status                    VARCHAR(16) NOT NULL DEFAULT 'active',
    owner_flag                TINYINT AS (CASE WHEN role = 'owner' THEN 1 END) VIRTUAL,

    UNIQUE (username),
    UNIQUE KEY uq_admin_users_email (email),
    UNIQUE KEY uq_admin_users_owner (owner_flag),
    CONSTRAINT chk_admin_users_role CHECK (role IN ('user', 'moderator', 'admin', 'owner')),
    CONSTRAINT chk_admin_users_status CHECK (status IN ('active', 'disabled'))
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
    key_prefix       CHAR(8) NULL,
    expires_at       TIMESTAMP NULL,

    CONSTRAINT fk_admin_api_keys_admin_user
        FOREIGN KEY (admin_user_id) REFERENCES admin_users(id) ON DELETE CASCADE,
    UNIQUE (key_hash),
    KEY idx_admin_api_keys_admin_user_id (admin_user_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- What each API key may do (API.9, v13): `read`, `generate`, `upload`,
-- `admin` (docs/design/api-design-standards.md section 5). A row per scope,
-- not a JSON column or a single scope column, so a key can hold
-- "upload and generate but not admin". No CHECK on `scope`: a new scope
-- needs no migration. `key_prefix` (the key's first eight characters) is
-- shown in the key list and the activity log so a key can be told apart
-- without the key; `expires_at` makes an expired key behave like a revoked one.
CREATE TABLE IF NOT EXISTS admin_api_key_scopes (
    key_id           BIGINT UNSIGNED NOT NULL,
    scope            VARCHAR(32) NOT NULL,

    PRIMARY KEY (key_id, scope),
    CONSTRAINT fk_admin_api_key_scopes_key
        FOREIGN KEY (key_id) REFERENCES admin_api_keys(id) ON DELETE CASCADE
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

-- v2 (SEC.1, SEC.21) kept failed-login counters here; v8 (SEC.30) moved
-- them to Redis (`planetgen/admin/throttle.py`'s RedisStore). The old
-- table is dropped, so a deployment loses the counts once.
DROP TABLE IF EXISTS login_throttle;

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

-- v5 (PERF.8): the generation work queue (`planetgen/queue/work.py`).
-- `work_jobs` is one row per run that queued work (a `planetgen`
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
    control          VARCHAR(16) NULL,        -- ADM.10: 'pause' or 'cancel', asked from the web

    KEY idx_work_jobs_state (state),
    KEY idx_work_jobs_created_at (created_at),
    KEY idx_work_jobs_root (root_id),
    CONSTRAINT fk_work_jobs_parent FOREIGN KEY (parent_id) REFERENCES work_jobs(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v7 (ADM.12): every job is a tree. A `planetgen` run is the root
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

-- v10 (DB.13): a job's command line without --mysql-* options, one row per
-- argument in order (it was one JSON text column, `work_jobs.argv`). A task's
-- `result` summary, written but never read, is gone.
CREATE TABLE IF NOT EXISTS work_job_args (
    job_id           VARCHAR(32) NOT NULL,
    position         INT NOT NULL,
    value            VARCHAR(1024) NOT NULL,
    PRIMARY KEY (job_id, position),
    CONSTRAINT fk_work_job_args_job FOREIGN KEY (job_id) REFERENCES work_jobs(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

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

-- v6 (PERF.3, PERF.10), v12 (PERF.32): how fast this server generates and
-- how much space a galaxy takes (`planetgen/generation/stats.py`).
-- `generation_stats` is one row per kind of task ("sector" fill, a
-- "scatter" layer of bright stars, "phenomena"), worker count and
-- log-scale density bucket (two per decade from 0.01, open-ended upward),
-- each a decaying average over every task that finished in it with that
-- many workers. `version_key` is the release, Python, OS and architecture
-- that measured it (`galaxy/version_key.py`): rows of any other key are
-- deleted by the next run, so a new release starts from nothing (Boss,
-- PERF.32), and an admin can delete them all. A benchmark's rows (PERF.31)
-- use a "bench:" kind prefix so they never feed a live estimate.
-- `generation_size` is one row per galaxy database: its bytes per star
-- system, measured from its own tables after each run. Bulk generation
-- reads both for its size and time estimate. A galaxy reset keeps them
-- (they describe the server).
CREATE TABLE IF NOT EXISTS generation_stats (
    kind                 VARCHAR(32) NOT NULL,     -- sector, scatter, phenomena, bench:...
    workers              INT NOT NULL DEFAULT 1,   -- worker processes of the runs it averages
    bucket               INT NOT NULL,             -- floor(2 * log10(density / 0.01))
    version_key          CHAR(22) NOT NULL DEFAULT '',
    density_low          DOUBLE NOT NULL,
    density_high         DOUBLE NOT NULL,
    samples              BIGINT UNSIGNED NOT NULL DEFAULT 0,
    seconds_per_task     DOUBLE NOT NULL DEFAULT 0,  -- worker wall time
    seconds_per_system   DOUBLE NOT NULL DEFAULT 0,
    systems_per_task     DOUBLE NOT NULL DEFAULT 0,
    stars_per_system     DOUBLE NOT NULL DEFAULT 0,
    max_density          DOUBLE NOT NULL DEFAULT 0,
    updated_at           DATETIME(6) NOT NULL,

    PRIMARY KEY (kind, workers, bucket)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS generation_size (
    database_name        VARCHAR(64) NOT NULL PRIMARY KEY,
    bytes_per_system     DOUBLE NOT NULL,
    systems              BIGINT UNSIGNED NOT NULL,
    total_bytes          BIGINT UNSIGNED NOT NULL,
    measured_at          DATETIME(6) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v9 (GEN.70): the naming key of each galaxy database, one row per
-- database. It is drawn from the galaxy seed when the galaxy is planned
-- and an admin can change it; the objects the phoneme codec names
-- (`planetgen/names/naming_key.py`) get their names from their ID and
-- this key, so changing it renames them all without rewriting a row.
-- `codec_version` is the codec's version when the key was set. A galaxy
-- reset keeps the row until a new plan draws a new key.
CREATE TABLE IF NOT EXISTS galaxy_naming (
    database_name        VARCHAR(64) NOT NULL PRIMARY KEY,
    naming_key           CHAR(8) NOT NULL,        -- 8 uppercase hex digits
    codec_version        SMALLINT UNSIGNED NOT NULL,
    drawn_at             DATETIME(6) NOT NULL,
    changed_at           DATETIME(6) NULL,        -- the last admin change
    changed_by           VARCHAR(64) NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v11 (OPS.13): the version key each update ran under, per galaxy database,
-- newest last. An update adds one row per galaxy (the galaxy's seed, the
-- key `planetgen.galaxy.version_key` computes, the release, and the SHA-256
-- of requirements.lock) and keeps only the last 10 per database
-- (`planetgen/galaxy/version_history.py`), so an admin can see which
-- release/environment a seed was last known to make the same galaxy under.
-- v18 (OPS.58): the row also holds the generator epoch, the battery digest
-- (and its cases) and the environment, and an update that finds all of them
-- as in the newest row adds none: it sets that row's `last_seen_at`.
CREATE TABLE IF NOT EXISTS version_key_history (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    database_name        VARCHAR(64) NOT NULL,
    galaxy_seed          CHAR(32) NOT NULL,       -- 32 uppercase hex digits
    version_key          CHAR(22) NOT NULL,       -- 22 uppercase hex digits
    planetgen_version    VARCHAR(32) NOT NULL,
    requirements_sha256  CHAR(64) NULL,           -- NULL: no requirements.lock beside the code
    recorded_at          DATETIME(6) NOT NULL,
    generator_epoch      SMALLINT UNSIGNED NULL,  -- v18 (OPS.58): NULL on a row recorded before it
    battery_digest       CHAR(64) NULL,           -- v18: the generator battery's overall digest (OPS.57)
    battery_cases        JSON NULL,               -- v18: {case: full digest}, to name the first case that differs
    environment          JSON NULL,               -- v18: libc, numpy, astropy, scipy and scikit-image versions
    last_seen_at         DATETIME(6) NULL,        -- v18: the last update that ran unchanged (NULL: recorded_at)
    INDEX idx_version_key_history_database (database_name, id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v14 (PERF.56): how long each stage of a generation run took, with the
-- settings it ran with (`planetgen/generation/stages.py` writes a row as a
-- stage ends). A stage a run skipped gets a row too, with `skipped` set and
-- its reason and no time. `settings` holds the arguments that shape the
-- stage (mass limit, luminosity floor, workers, radius ...) and `metrics`
-- what it did (layers visited, layers that placed something, objects), both
-- as JSON text, so a later estimate can pick the rows run with the same
-- settings. Rows of another `version_key` are not used for estimates.
CREATE TABLE IF NOT EXISTS generation_stage_runs (
    id                   BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    database_name        VARCHAR(64) NOT NULL,
    command              VARCHAR(16) NOT NULL,     -- galaxy, plan
    stage_key            VARCHAR(24) NOT NULL,
    stage_n              SMALLINT UNSIGNED NOT NULL,
    stage_total          SMALLINT UNSIGNED NOT NULL,
    label                VARCHAR(80) NOT NULL,
    skipped              TINYINT(1) NOT NULL DEFAULT 0,
    skip_reason          VARCHAR(200) NULL,
    started_at           DATETIME(6) NOT NULL,
    finished_at          DATETIME(6) NOT NULL,
    seconds              DOUBLE NOT NULL DEFAULT 0,
    workers              INT NOT NULL DEFAULT 1,
    settings             TEXT NOT NULL,
    metrics              TEXT NOT NULL,
    version_key          CHAR(22) NOT NULL DEFAULT '',
    INDEX idx_generation_stage_runs_key (stage_key, skipped, id),
    INDEX idx_generation_stage_runs_database (database_name, id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v15 (DB.58): one-time tokens for the account flows (USR.28: invite, first
-- password, reset, email change and the two links of an owner transfer). Only
-- the SHA-256 hex digest of the token is stored. `expires_at` ends it, `used_at`
-- is set once when it is redeemed (the redeem is one UPDATE that requires
-- `used_at IS NULL`, so a token works once), and `meta` holds what the purpose
-- needs (an invite's use count, note and invitee email). The tokens of an
-- account go with the account.
CREATE TABLE IF NOT EXISTS account_tokens (
    id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    account_id    BIGINT UNSIGNED NOT NULL,
    purpose       VARCHAR(24) NOT NULL,
    token_hash    CHAR(64) NOT NULL,
    created_at    TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    expires_at    TIMESTAMP NULL,
    used_at       TIMESTAMP NULL,
    meta          JSON NULL,

    UNIQUE KEY uq_account_tokens_hash (token_hash),
    KEY idx_account_tokens_account (account_id, purpose),
    CONSTRAINT fk_account_tokens_account
        FOREIGN KEY (account_id) REFERENCES admin_users(id) ON DELETE CASCADE,
    CONSTRAINT chk_account_tokens_purpose CHECK (purpose IN
        ('invite', 'first_password', 'reset', 'email_change', 'owner_transfer_old', 'owner_transfer_new'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v15 (DB.58): who redeemed an invite link. A link with several uses has one
-- row per account that registered through it. Removed with the invite token
-- or with the account.
CREATE TABLE IF NOT EXISTS invite_uses (
    token_id      BIGINT UNSIGNED NOT NULL,
    account_id    BIGINT UNSIGNED NOT NULL,
    used_at       TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    PRIMARY KEY (token_id, account_id),
    KEY idx_invite_uses_account (account_id),
    CONSTRAINT fk_invite_uses_token
        FOREIGN KEY (token_id) REFERENCES account_tokens(id) ON DELETE CASCADE,
    CONSTRAINT fk_invite_uses_account
        FOREIGN KEY (account_id) REFERENCES admin_users(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v16 (DB.62): bookkeeping for the uploads of generated sectors from a remote
-- run (API.11; docs/design/api-design-standards.md section 7.3). Nothing
-- reads or writes these tables until API.11. One `upload_runs` row is one
-- remote run: its owner (`owner_id`, an account) and the API key it used
-- (`key_id`, NULL once the key is removed), the galaxy database it uploads
-- to, its `state` (set by API.11) and when it began and was last heard from.
CREATE TABLE IF NOT EXISTS upload_runs (
    id                BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    owner_id          BIGINT UNSIGNED NOT NULL,
    key_id            BIGINT UNSIGNED NULL,
    galaxy_database   VARCHAR(64) NOT NULL,
    state             VARCHAR(16) NOT NULL DEFAULT 'open',
    created_at        TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    last_contact_at   TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),

    KEY idx_upload_runs_owner (owner_id, state),
    KEY idx_upload_runs_galaxy (galaxy_database, state),
    CONSTRAINT fk_upload_runs_owner
        FOREIGN KEY (owner_id) REFERENCES admin_users(id) ON DELETE CASCADE,
    CONSTRAINT fk_upload_runs_key
        FOREIGN KEY (key_id) REFERENCES admin_api_keys(id) ON DELETE SET NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v16 (DB.62): the gzip batches a run has sent, one row per (run, sequence
-- number): the SHA-256 of the received bytes and their size, so a resent batch
-- with the same hash is answered from the record and a different one is refused.
CREATE TABLE IF NOT EXISTS upload_batches (
    run_id        BIGINT UNSIGNED NOT NULL,
    seq           INT UNSIGNED NOT NULL,
    sha256        CHAR(64) NOT NULL,
    bytes         BIGINT UNSIGNED NOT NULL,
    received_at   TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),

    PRIMARY KEY (run_id, seq),
    CONSTRAINT fk_upload_batches_run
        FOREIGN KEY (run_id) REFERENCES upload_runs(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v16 (DB.62): the units (a sector with what is in it) of a run, each written to
-- the galaxy only when complete: `open` while batches arrive, `sealed` when the
-- client says it is whole, `verified` when the counts and hash match the
-- manifest, then `finalized` (saved) or `rejected` with its reasons in `errors`.
CREATE TABLE IF NOT EXISTS upload_units (
    run_id          BIGINT UNSIGNED NOT NULL,
    unit_key        VARCHAR(64) NOT NULL,
    state           VARCHAR(12) NOT NULL DEFAULT 'open',
    manifest_hash   CHAR(64) NULL,
    errors          JSON NULL,
    updated_at      TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3) ON UPDATE CURRENT_TIMESTAMP(3),

    PRIMARY KEY (run_id, unit_key),
    KEY idx_upload_units_state (state),
    CONSTRAINT fk_upload_units_run
        FOREIGN KEY (run_id) REFERENCES upload_runs(id) ON DELETE CASCADE,
    CONSTRAINT chk_upload_units_state CHECK (state IN ('open', 'sealed', 'verified', 'finalized', 'rejected'))
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v16 (DB.62): which run holds a sector's address. One claim per sector address
-- in a galaxy database (`galaxy_database` repeats the run's, so the key can say
-- so), so a second run asking for a claimed sector is refused.
CREATE TABLE IF NOT EXISTS upload_claims (
    galaxy_database   VARCHAR(64) NOT NULL,
    shell_index       INT NOT NULL,
    band_index        INT NOT NULL,
    wedge_index       INT NOT NULL,
    run_id            BIGINT UNSIGNED NOT NULL,
    state             VARCHAR(16) NOT NULL DEFAULT 'claimed',
    claimed_at        TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),

    PRIMARY KEY (galaxy_database, shell_index, band_index, wedge_index),
    KEY idx_upload_claims_run (run_id),
    CONSTRAINT fk_upload_claims_run
        FOREIGN KEY (run_id) REFERENCES upload_runs(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v16 (DB.62): the row ids a run may use, from the galaxy's id blocks (API.10):
-- one contiguous range per table and lease, `first_id` to `last_id` inclusive.
CREATE TABLE IF NOT EXISTS upload_id_ranges (
    run_id       BIGINT UNSIGNED NOT NULL,
    table_name   VARCHAR(64) NOT NULL,
    first_id     BIGINT UNSIGNED NOT NULL,
    last_id      BIGINT UNSIGNED NOT NULL,

    PRIMARY KEY (run_id, table_name, first_id),
    CONSTRAINT fk_upload_id_ranges_run
        FOREIGN KEY (run_id) REFERENCES upload_runs(id) ON DELETE CASCADE,
    CONSTRAINT chk_upload_id_ranges_order CHECK (last_id >= first_id)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- v17 (ADM.125): Generate recipes, a named and saved Generate form kept per
-- admin. A recipe's form is the same draft "Run again" keeps for a job, one
-- `generate_recipe_fields` row per `[name, value]` pair in the order the form
-- sent them (its `action` field says which task). An admin's names are
-- unique, and the recipes go with the account and the fields with the recipe.
CREATE TABLE IF NOT EXISTS generate_recipes (
    id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    admin_id      BIGINT UNSIGNED NOT NULL,
    name          VARCHAR(80) NOT NULL,
    created_at    TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at    TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,

    UNIQUE KEY uq_generate_recipes_name (admin_id, name),
    CONSTRAINT fk_generate_recipes_admin
        FOREIGN KEY (admin_id) REFERENCES admin_users(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

CREATE TABLE IF NOT EXISTS generate_recipe_fields (
    recipe_id     BIGINT UNSIGNED NOT NULL,
    position      SMALLINT UNSIGNED NOT NULL,
    field_name    VARCHAR(100) NOT NULL,
    field_value   VARCHAR(4000) NOT NULL,

    PRIMARY KEY (recipe_id, position),
    CONSTRAINT fk_generate_recipe_fields_recipe
        FOREIGN KEY (recipe_id) REFERENCES generate_recipes(id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
