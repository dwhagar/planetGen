### Changed
- Generation runs on the shell grid: regions and spans are shells, bands and wedges (`galaxy/region.py`, `span.py`); the bright-star and phenomenon scatters draw one shell at a time from a stream of their own, each band with its own density bound (`object_first.band_majorants`); the nebula field's extent comes from the outline (`nebula_field.field_extent_pc`).
- `planetgen galaxy` takes `--shell I` (with `--band J` and `--wedge K` for one band or one sector), spans `--shells/--bands/--wedges`, and random starts `--max-shell`; `plan` takes `--max-shell`; `fingerprint` and `check-db` take `--shell`. `limits.MAX_GENERATE_SHELL` replaces `MAX_GENERATE_RING`. The cylinder fill (`--cylinder-sectors`, `--cylinder-layers`), `--column`, `--block-layer`, `--ring`, `--layer`, `--slot`, `--rings`, `--layers`, `--slots` and `--max-ring` are gone; `--block` is refused until the Galaxy Map's ladder is rebuilt on shells (MAP.205).
- The neighborhood fill is a ball: ADM.30's "cylinder height = X layers" no longer applies.

TODO items resolved: GEN.227 (in part: scatters, regions, spans, flags, limits) / Added: none
