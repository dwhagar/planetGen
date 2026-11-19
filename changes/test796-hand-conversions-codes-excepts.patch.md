**Headline:** Four hand mass conversions use the unit registry, no TODO code is shown on the property sheet or a limit's text, and two quiet handlers log at debug.

### Fixed
- `generation/cosmic_dust.py`, `physics/disk_truncation.py`, `physics/envelope_capture.py` and `physics/influx.py` multiplied or divided by `SOLAR_MASS_TO_KG` or `EARTH_MASS_TO_KG` by hand (`test_mass_conversions`); they use `unit_registry.converter` now. The converter applies the same single multiply or divide, so the numbers are the same bit for bit (checked on 600,000 random masses).
- The property sheet (`describe/sheet.py`), the class fit (`physics/class_fit.py`) and the protection limits (`physics/protection_limits.py`) showed item codes (DB.133, GEN.1752, GEN.1231, GEN.1737, GEN.1772) in strings a reader sees; the codes are gone and the sentences are unchanged otherwise (`test_no_todo_codes_shown`).
- `util/detmath._kernel` and `web/system_page.generate_in_process` swallowed an exception without a word; each now says it at debug level (`log.swallowed`), behaviour unchanged (`test_broad_except_logs`).

generation-output: unchanged
