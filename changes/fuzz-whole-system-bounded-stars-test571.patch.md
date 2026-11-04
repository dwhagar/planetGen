**Headline:** The whole-system fuzz test no longer draws the slow giant and supergiant star classes and runs 20 examples, so it finishes instead of timing out.

### Changed

- `test_random_valid_configs_generate_sane_systems` ran past 2900 s on main: a star drawn from the Yerkes classes IA+, IA, IAB, IB, II or III costs 10 to 60 s of per-body air and climate work, and half the examples name a star type. It now runs 20 examples instead of 80 (still scaled up in the deep profile), and its configs draw a star type from the other classes (420 of 840 types). Every spectral letter and subclass of the heavy classes stays covered by `test_every_star_type_generates_a_sane_system`, and every sanity check of the whole-system test is unchanged. The test also stopped asserting that `HABITABLE_WORLD=False` leaves no life-class body (the bar is only the request's conditions in the ecosphere, GEN.916); it asserts the system makes no habitable request. Tests only.

generation-output: unchanged
