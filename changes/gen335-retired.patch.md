**Headline:** The TODO list retires the single definition of the stored point columns and the galaxy-grid fuzz bug.

### Changed
- **GEN.335 and TEST.192 are done (#1883).** The stored point of an object is defined once, with a test that lists every object table, and the fuzz test's tolerance now allows for the rounding of the halo subtraction. The migration order is now DB.55 0101, DB.54 0102, GEN.336 0103, DB.53 0104; belts and asteroid fields store a stated zero mass until a mass model is asked for, which NAV.62 notes.

TODO items resolved: GEN.335, TEST.192 / Added: none
