**Headline:** The admin's one-off system page takes its random seed through the draw module again, which fixes the red `test_generation_draws_only_through_the_draw_module`.

OPS.244 slice A moved the one-off system run into the supervisor's warm process and, on the way, drew its seed with `secrets.randbits` in `generation/run_system.py` and `web/system_page.py`. Generation code may not import `secrets`. `util/draw.py`, the one module that is allowed to, now has `fresh_seed()` for a draw that is not part of any run, and both callers use it.

generation-output: unchanged
generation-affects: none
