"""The galaxy: its sector grid (`spherical.py`: shells, bands and wedges, with
the designations and the core ball), density model (`density.py`) and outline
(`skeleton.py`), the seed and version key (`seed.py`, `version_key.py`,
`version_history.py`, `version_check.py`), object IDs (`uid.py`,
`object_uid.py`, `objectref.py`), sectors (`sector.py`), the molecular cloud
field and nebula shapes (`nebula_field.py`, `nebula_shape.py`), positions and
motion (`frame.py`, `galactic_orbit.py`, `system_position.py`) and navigation
(`navigation.py`, `nav_graph.py`, `straight_line.py`, `keepout.py`).

The Galaxy Map's block ladder on shells is `drill.py` (blocks of 243 to 3 shells,
docs/design/spherical-block-ladder.md), its data windows `window.py`, its cube
tiles `viewport.py` and the regions the generation commands fill `region.py`.
`geometry.py` holds only what is left of the old cylinder grid: the offset of a
binary companion."""
