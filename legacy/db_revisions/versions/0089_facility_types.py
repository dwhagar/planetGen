"""Facility types, their typed fields and values, and the new facility columns (DB.28, POP.10).

Schema v89. Creates `facility_types`, `facility_type_fields` and `facility_field_values`, and adds `type_id`,
`affiliation_polity_id`, `affiliation_name`, the five ratings (`rating_crime`, `rating_housing`,
`rating_resources`, `rating_maintenance`, `rating_health`, 1 Green to 3 Red), `ratings_source`, `ratings_note`
and `ratings_set_at` to `facilities`. All NULL for facilities saved before v89. Each step is skipped when it is
already done, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0089"
down_revision = "0088"
branch_labels = None
depends_on = None

_TABLE_OPTIONS = "ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"
_RATINGS = ("rating_crime", "rating_housing", "rating_resources", "rating_maintenance", "rating_health")


def _has_table(connection, table):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name"), {"table_name": table}).scalar()


def _has_column(connection, table, column):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar()


def upgrade():
    connection = op.get_bind()
    if not _has_table(connection, "facility_types"):
        op.execute(f"""
            CREATE TABLE facility_types (
                id           BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
                name         VARCHAR(255) NOT NULL,
                kind         VARCHAR(16) NOT NULL,
                icon         VARCHAR(64),
                description  TEXT,
                sort_order   INT NOT NULL DEFAULT 0,
                active       TINYINT(1) NOT NULL DEFAULT 1,
                UNIQUE KEY uq_facility_types_name (name),
                KEY idx_facility_types_kind (kind),
                CONSTRAINT chk_facility_types_kind
                    CHECK (kind IN ('colony', 'outpost', 'mining-colony', 'station', 'starbase')),
                CONSTRAINT chk_facility_types_active CHECK (active IN (0, 1))
            ) {_TABLE_OPTIONS}""")
    if not _has_table(connection, "facility_type_fields"):
        op.execute(f"""
            CREATE TABLE facility_type_fields (
                id            BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
                type_id       BIGINT UNSIGNED NOT NULL,
                `key`         VARCHAR(32) NOT NULL,
                label         VARCHAR(255) NOT NULL,
                datatype      VARCHAR(8) NOT NULL,
                unit          VARCHAR(32),
                required      TINYINT(1) NOT NULL DEFAULT 0,
                min_value     DOUBLE,
                max_value     DOUBLE,
                choices       TEXT,
                default_text  VARCHAR(255),
                sort_order    INT NOT NULL DEFAULT 0,
                UNIQUE KEY uq_facility_type_fields_key (type_id, `key`),
                CONSTRAINT chk_facility_type_fields_datatype
                    CHECK (datatype IN ('int', 'decimal', 'text', 'bool', 'choice')),
                CONSTRAINT chk_facility_type_fields_required CHECK (required IN (0, 1)),
                CONSTRAINT fk_facility_type_fields_type
                    FOREIGN KEY (type_id) REFERENCES facility_types(id) ON DELETE CASCADE
            ) {_TABLE_OPTIONS}""")
    if not _has_column(connection, "facilities", "type_id"):
        columns = ", ".join(
            ["ADD COLUMN type_id BIGINT UNSIGNED NULL AFTER uid",
             "ADD COLUMN affiliation_polity_id BIGINT UNSIGNED NULL AFTER type_id",
             "ADD COLUMN affiliation_name VARCHAR(255) NULL AFTER affiliation_polity_id"]
            + [f"ADD COLUMN {rating} TINYINT UNSIGNED NULL AFTER {after}" for rating, after in zip(
                _RATINGS, ("affiliation_name",) + _RATINGS)]
            + ["ADD COLUMN ratings_source VARCHAR(16) NULL AFTER rating_health",
               "ADD COLUMN ratings_note TEXT NULL AFTER ratings_source",
               "ADD COLUMN ratings_set_at TIMESTAMP NULL AFTER ratings_note"])
        constraints = ", ".join(
            [f"ADD CONSTRAINT chk_facilities_{rating} CHECK ({rating} BETWEEN 1 AND 3)" for rating in _RATINGS]
            + ["ADD KEY idx_facilities_type (type_id)", "ADD KEY idx_facilities_affiliation (affiliation_polity_id)"]
            + [f"ADD KEY idx_facilities_{rating} ({rating})" for rating in _RATINGS]
            + ["ADD CONSTRAINT fk_facilities_type FOREIGN KEY (type_id) REFERENCES facility_types(id)"
               " ON DELETE RESTRICT",
               "ADD CONSTRAINT fk_facilities_affiliation_polity FOREIGN KEY (affiliation_polity_id)"
               " REFERENCES polities(id) ON DELETE SET NULL"])
        op.execute(f"ALTER TABLE facilities {columns}")
        op.execute(f"ALTER TABLE facilities {constraints}")
    if not _has_table(connection, "facility_field_values"):
        op.execute(f"""
            CREATE TABLE facility_field_values (
                facility_id  BIGINT UNSIGNED NOT NULL,
                field_id     BIGINT UNSIGNED NOT NULL,
                value_int    BIGINT,
                value_dec    DOUBLE,
                value_text   VARCHAR(255),
                value_bool   TINYINT(1),
                PRIMARY KEY (facility_id, field_id),
                KEY idx_facility_field_values_field (field_id),
                CONSTRAINT fk_facility_field_values_facility
                    FOREIGN KEY (facility_id) REFERENCES facilities(id) ON DELETE CASCADE,
                CONSTRAINT fk_facility_field_values_field
                    FOREIGN KEY (field_id) REFERENCES facility_type_fields(id) ON DELETE RESTRICT
            ) {_TABLE_OPTIONS}""")
