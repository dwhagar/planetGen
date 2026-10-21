**Headline:** The test that each class draws from a stream of its own now scatters a hundredth of the rows, so it no longer takes ten gigabytes and gets a gate worker killed.

### Fixed
- Tests: `test_each_class_draws_from_a_stream_of_its_own_so_the_black_holes_move_no_other_row` listed ten million scatter rows twice (10.5 GB peak, 164 s), and in a four-worker gate the box ran out of memory and killed a worker (TEST.256). It now scatters 40 expected per cell (100 000 rows, 0.24 GB, 2 s) and still proves the same two things: no other class's row moves when the black hole rate triples, and the black holes grow by more than half.

TODO items resolved: TEST.256 / Added: none
