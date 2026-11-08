**Headline:** A request for a planet now names conditions (a kind, an air, a range for any quantity the classes read); a class letter is only a short way to write a class's conditions, and the finished planet is judged on them, never drawn toward a letter (GEN.688).

### Changed

- New `generation/condition_request.py`: `from_letter` (a class's own conditions from `planet_class_fit.json`; a letter two kinds share, T, is either), `from_conditions` (kind, air, `[low, high]` ranges for the quantities the classes read), `meets` (a finished planet or moon meets every condition of an option, with hard edges) and `wants_ecosphere` (a rocky body whose temperature must lie between freezing and 400 K is placed in the ecosphere; stand-in CRQ1, section 6.218).
- Slots take `"conditions"` beside `"planet_class"` (never both); `SystemConfig.FORCE_CONDITIONS` joins `FORCE_CLASS`. `StarSystem.build_planet` draws the generator's own planets (up to 40 at a place) until one meets the request; the log-uniform mass step (`CLASS_REQUEST_MASS_EARTH`, `CLASS_REQUEST_FREE_ATTEMPTS`, registry PA14), the request-only `late_core` branch (`CLASS_REQUEST_LATE_CORE_EVERY`, registry EC7; the free generator never takes it) and the slot distance snap by the life-bearing class list are removed.
- A request no longer draws the whole system again: the class part of `planets_regenerated` and `planets_regenerated_after_age` is gone, and the system reports the request in `unmet_requirements` when no planet meets it. The rare classes (P, M, H, O, R, and V without `late_core`) stay requestable and come out "unmet" where the free draw cannot make them (Boss, GEN.1333).
- `FORCE_CLASS` is judged by the class's conditions, not by the letter the classifier gives the planet.

generation-output: unchanged
