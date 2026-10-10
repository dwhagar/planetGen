### Added
- The object ID layout (GEN.170): `galaxy/object_uid.py` packs, unpacks, prints and parses the 80-bit ID of birth sector, serial and body number (20 hex digits, `BINARY(10)`), and picks a longer layout (96 or 128 bits) for a galaxy whose bounds do not fit. Nothing uses it yet; the schema (DB.20) and the fill (GEN.171) come next.
