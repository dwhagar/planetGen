**Headline:** the job engine design note now describes the engine as it is on main, with RQ gone, and says what is proven and what is not.

### Changed
- OPS.252 (the rest): `docs/design/operations/job-engine.md` is brought up to date for what merged since its first slice: the stage table (stages 3 to 6: OPS.235 RQ removed, OPS.239, 240, 244 done, OPS.249, 250, 251, 261, 271 done), section 3 (failures reach the engine as they are), section 5 (sectors one at a time with a spawn row and a window, claims, time limits and the split follow the profile of a proven kind), section 6 (the controller's edges are learned), section 7 ("The old queue is gone": what stays under the old names and is not a queue), section 8 (the cold one-off system request measured on a quiet machine: 52.6 s first, 19.5 to 28.4 s after), section 10 (what is proven: OPS.225 slices 1 to 3, OPS.274; what is not) and the code map.

generation-output: unchanged
generation-affects: none
