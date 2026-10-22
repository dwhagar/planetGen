**Headline:** Every persisted field of the stellar classes now declares the unit it is held in, and `star.quantity("mass")` answers in any other.

### Added
- **Field declarations** (GEN.570, first slice): `MassObject.FIELD_UNITS` maps each persisted attribute to a registry unit
  (or a plain word such as `text`, `flag`, `ratio`), `SERIALIZABLE_FIELDS` is its names, and `obj.quantity(name)` builds
  the `Quantity` in the declared native unit. Declared for `Star`, `BlackHole`, `NeutronStar` and `BinaryStarProxy`; the
  shared spin, activity, disk and galactic-motion groups are declared once in their modules. The instance still holds a
  plain float, so generation does no more work than before.
- The unit registry holds `1/s`, `1/yr` and `1/Gy` as frequency units; the column registry's local copies are gone.
- `tests/test_field_units.py` fails on an undeclared persisted field, a name whose suffix disagrees with its declaration,
  or a declared unit that measures something other than its column.
