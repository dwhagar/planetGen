**Headline:** `SectorCell.contains` now counts a point on a far cap's colatitude face and on the last wedge's end face as inside.

### Fixed
- **Two points `SectorCell.sample` draws exactly on a face were taken for outside the cell (found by the deep fuzz profile).** On a cap at a far shell (shell 1,000,000, a sector 1e6 pc out) the colatitude of a point is only as good as its coordinates, which are rounded at a few 1e-16 of the radius, over the point's distance from the axis: 2.4e-11 rad off the face, where the check allowed 1e-12. The allowance now adds 1e-15 over the sine of the cell's colatitude. And a point on the end face of the last wedge of a band, at longitude 2 pi, reads back as 6e-16 and was taken for outside the wedge. Both are caught by `contains` only; sampling and storage are unchanged, so no generated sector moves. `test_the_end_face_of_a_last_wedge_is_inside_though_its_longitude_reads_back_as_zero` and an `@example` on `test_cell_samples_on_every_face_are_inside` fail without the change.

TODO items resolved: none. Added: none.
