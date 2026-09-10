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
--     `check_password_hash`. See `stellarObjects/adminAuth.py`.
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
-- that performs the action (`stellarObjects.adminAuth.record_audit`), not
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
