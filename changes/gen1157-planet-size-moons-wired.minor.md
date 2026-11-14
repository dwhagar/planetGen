**Headline:** A few giant planets now hold a moon as big as a planet, taken by capture, while the moons of a giant's own disk stay under the most a disk can make.

### Changed

- GEN.1157: a giant's moons are held to the most a circumplanetary disk can hold of its parent between them, 6.4e-4 of its mass (`planet_moons.disk_moon_ceiling_kg`, now read by `planet_kind.moon_mass_ceiling_kg`; the flat 1e-3 cap, `tuning.MOON_GIANT_PARENT_MAX_MASS_RATIO`, is removed). A Jupiter's moons now total at most 0.2 Earth masses; a super-Jovian of ten Jupiter masses can still hold a moon of a couple of Earth masses.
- About one giant in 70 also holds one moon of planet size taken by capture (a tenth of an Earth mass up to 1e-2 of the giant), drawn after its disk moons by `planet_moons.draw_exchange_moon` and placed past them. It has the new moon origin `exchange`: its orbit is a captured one (isotropic outside the Kozai band, so about half are retrograde, with a thermal eccentricity), and its body is made as any planet's is, so its class, air and life follow from its own conditions.
- Stand-ins PSM1 to PSM4 (section 6.276) are read now; CC2 names the new cap. No migration: `moon_origin` is VARCHAR(8) and `exchange` is eight characters.
- Where a giant has no room past its last moon, or the drawn body is too dense or light for its parent to hold, the capture moon is skipped.
