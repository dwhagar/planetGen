### Fixed

- Sector and galaxy runs can make a feature more or less common instead of forcing it on every system: `--prevalence comets=+50` gives 1.5 times as many systems with comets, `--prevalence habitable_world=-100` none. It works for every former forcing option (habitable world, asteroid belt, comets, large star, moons, max planets, intelligent life, binary, wide binary, planets), and each system stores the setting with its config (database schema v54).
- `-large_star` on a single system now really keeps out a large star; it used to make no difference.
