**Headline:** A rocky planet under a deep envelope whose haze keeps its 1 bar level colder than its albedo alone gives now stores a temperature that matches its stored albedo and greenhouse quantity.

### Changed
- **The balance of a hazy envelope is kept in the albedo (Train 28 sweep).** The solver found a negative greenhouse quantity for a body whose haze takes more sunlight from the ground than its greenhouse returns (the anti-greenhouse of Titan); the stored quantity is never negative, and `air_final` clipped it to 0.0, which left the stored temperature up to 0.26 percent off the formula on the stored albedo and greenhouse (one body of the sweep: 180.77 K against 181.25 K). The greenhouse quantity stays 0.0 and the Bond albedo carries the balance (`climate_draws.albedo_for_temperature`: `1 - A = 4 sigma T^4 / S`), so the temperature does not move and the pair reproduces it. Only a giant or envelope body whose temperature is below its airless one changes (its albedo rises, 0.3383 to 0.3453 for the body found); nothing else does.

generation-output: changed
generation-affects: fill
