"""The orbital event record (DB.56).

Schema v118. Two new tables: `orbital_events`, one row per orbital event (a collision outcome, merger, split, ejection,
capture, unbind or warning) with the numbers it was decided from, the conservation check, the products and the flags of
docs/design/collisions-and-mergers.md section 7, and `orbital_event_bodies`, the bodies that went into it and came out of
it, one row each. The revision only makes the tables: nothing writes them yet (GEN.197's `record_event` does) and no
view reads them (GEN.143). A step already done is left alone, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0118"
down_revision = "0117"
branch_labels = None
depends_on = None

TABLES = {
    "orbital_events": """
        CREATE TABLE IF NOT EXISTS orbital_events (
            id                          BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
            run_id                      BIGINT UNSIGNED NULL,
            kind                        VARCHAR(16) NOT NULL,
            outcome_code                VARCHAR(24) NULL,
            regime                      VARCHAR(24) NULL,
            event_time_unix             DOUBLE NOT NULL,
            detected_time_unix          DOUBLE NOT NULL,
            seed_string                 VARCHAR(255) NOT NULL,
            separation_start_km         DOUBLE NULL,
            separation_impact_km        DOUBLE NULL,
            impact_parameter_km         DOUBLE NULL,
            critical_impact_parameter_km DOUBLE NULL,
            v_inf_kms                   DOUBLE NULL,
            v_esc_kms                   DOUBLE NULL,
            v_imp_kms                   DOUBLE NULL,
            q_r_j_per_kg                DOUBLE NULL,
            q_star_rd_j_per_kg          DOUBLE NULL,
            erosion_x                   DOUBLE NULL,
            interacting_mass_fraction   DOUBLE NULL,
            catastrophic_velocity_kms   DOUBLE NULL,
            window_s                    DOUBLE NULL,
            window_capped               TINYINT NOT NULL DEFAULT 0,
            focusing_radius_used        TINYINT NOT NULL DEFAULT 0,
            mass_error                  DOUBLE NULL,
            momentum_error              DOUBLE NULL,
            conservation_tolerance      DOUBLE NULL,
            impact_energy_j             DOUBLE NULL,
            heat_share                  DOUBLE NULL,
            lost_mass_kg                DOUBLE NULL,
            merged_mass_kg              DOUBLE NULL,
            new_radius_km               DOUBLE NULL,
            spin_breakup_ratio          DOUBLE NULL,
            deuterium_limit_mj          DOUBLE NULL,
            hydrogen_limit_sol          DOUBLE NULL,
            borderline                  TINYINT NOT NULL DEFAULT 0,
            containing_uid              BINARY(12) NULL,
            influence_body_count        INT UNSIGNED NULL,
            sectors_dirtied             INT UNSIGNED NULL,
            suppressed_duplicates       INT UNSIGNED NOT NULL DEFAULT 0,
            already_overlapping         TINYINT NOT NULL DEFAULT 0,
            numeric_warning             TINYINT NOT NULL DEFAULT 0,
            flags                       VARCHAR(255) NOT NULL DEFAULT '',
            recorded_at                 TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
            KEY idx_orbital_events_time (event_time_unix, id),
            KEY idx_orbital_events_kind (kind, event_time_unix, id),
            KEY idx_orbital_events_run (run_id, id),
            CONSTRAINT chk_orbital_events_kind CHECK (kind IN
                ('collision', 'merger', 'split', 'ejection', 'capture', 'unbind', 'warning')),
            CONSTRAINT chk_orbital_events_window_capped CHECK (window_capped IN (0, 1)),
            CONSTRAINT chk_orbital_events_focusing CHECK (focusing_radius_used IN (0, 1)),
            CONSTRAINT chk_orbital_events_borderline CHECK (borderline IN (0, 1)),
            CONSTRAINT chk_orbital_events_overlapping CHECK (already_overlapping IN (0, 1)),
            CONSTRAINT chk_orbital_events_numeric_warning CHECK (numeric_warning IN (0, 1))
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;""",
    "orbital_event_bodies": """
        CREATE TABLE IF NOT EXISTS orbital_event_bodies (
            event_id                    BIGINT UNSIGNED NOT NULL,
            role                        VARCHAR(12) NOT NULL,
            ordinal                     SMALLINT UNSIGNED NOT NULL,
            uid                         BINARY(12) NOT NULL,
            object_kind                 VARCHAR(24) NOT NULL,
            fate                        VARCHAR(12) NULL,
            name                        VARCHAR(255) NULL,
            object_class                VARCHAR(64) NULL,
            mass_kg                     DOUBLE NULL,
            radius_km                   DOUBLE NULL,
            composition                 VARCHAR(64) NULL,
            sector_id                   BIGINT UNSIGNED NULL,
            x_pc                        DOUBLE NULL,
            y_pc                        DOUBLE NULL,
            z_pc                        DOUBLE NULL,
            vx_kms                      DOUBLE NULL,
            vy_kms                      DOUBLE NULL,
            vz_kms                      DOUBLE NULL,
            admin_edited                TINYINT NOT NULL DEFAULT 0,
            note                        VARCHAR(255) NULL,
            PRIMARY KEY (event_id, role, ordinal),
            KEY idx_orbital_event_bodies_uid (uid, event_id),
            CONSTRAINT fk_orbital_event_bodies_event FOREIGN KEY (event_id) REFERENCES orbital_events(id) ON DELETE CASCADE,
            CONSTRAINT chk_orbital_event_bodies_role CHECK (role IN ('participant', 'product')),
            CONSTRAINT chk_orbital_event_bodies_fate CHECK (fate IS NULL OR fate IN ('destroyed', 'survived', 'absorbed', 'created')),
            CONSTRAINT chk_orbital_event_bodies_edited CHECK (admin_edited IN (0, 1))
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;""",
}


def upgrade():
    connection = op.get_bind()
    for statement in TABLES.values():    # `schema.sql` is applied first, so the tables are usually there already
        connection.execute(sa.text(statement))
