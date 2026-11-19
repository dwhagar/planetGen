**Headline:** An asteroid field's page now shows the field in 3D: its own bodies as a cloud at the places its age and speed put them, with the real body count stated under it.

### Added

- GEN.361: the view of an asteroid field (`web/maps/phenomenonrender.py` `field_view`, drawn by `static/phenomenonrender.js`, and the Sensor registry's entry for the kind). It draws 480 bodies of the field (the heaviest 12 larger), each at `field_dispersal.member_position_m` of the field's Object ID at its age and rms speed (a ball for a collision, a thin ring for a belt a supernova unbound), so the same field is the same cloud after a database rebuild. The line under it states how many bodies larger than 1 km the field really holds (from its mass and the size law of its makeup), that each point is drawn far larger than any of them, its column density end on (`column_density_kg_m2`), how long after it began to spread its 1 km pieces stopped colliding (`freeze_out_time_s`) and how far galactic shear has drawn its far side along its orbit (`shear_displacement_m`). A field older than the 26 Myr the shear is followed for says so.
- Four functions of `field_dispersal` that nothing read are now read by the page (the wiring register lines are deleted).

generation-output: unchanged
generation-affects: none
