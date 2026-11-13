**Headline:** An ocean world now needs 95 percent or more of its surface under liquid water (was 80 percent), and a temperate world runs up to 95 percent.

### Changed

- GEN.1570: the class edges for liquid cover move (Boss's edge, registry CL6, tier T6, L0): class O takes `ocean_fraction` 0.95 or more (was 0.8) and class M takes under 0.95 (was under 0.8); H is unchanged (under 0.2). The cover is the liquid sea since GEN.1567. A body at exactly 0.95 scores 1.0 in both and M takes it (the narrower description wins); from 0.951 up it is O, and a full ocean is O at 1.0 and M at 0. The fudge table was rebuilt from the new edges. No body in the census sits between 71 and 95 percent cover, so no stored letter is expected to move, but a body there now reads M where it read O.

generation-output: changed
generation-affects: fill
