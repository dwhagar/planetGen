# Example System Files

Run any of these with:

```bash
python planetGen.py --system-file examples/<file>.json
```

See [`JSON.md`](JSON.md) for the full format reference.

| File | System | Notes |
|---|---|---|
| `solar_system.json` | Our own solar system | Sun (G2V) through an asteroid belt to two gas giants and two ice giants. |
| `tatooine_system.json` | Tatooine (*Star Wars*) | A twin-sun (binary) system; Tatooine itself is class `H` (desert, minimal water), moonless. |
| `vulcan_system.json` | Vulcan (*Star Trek*) | Orbits a K-type star (40 Eridani); Vulcan is class `V` (super-Earth, high gravity) with one moon. |
| `krypton_system.json` | Krypton (*Superman*) | Orbits an aging red giant (Rao); Krypton is class `V` (dense, high-gravity world), `age: "old"`. |
| `arrakis_system.json` | Arrakis / Dune (*Dune*) | Orbits Canopus (`F0IB`, matching the real star's supergiant luminosity class); Arrakis is class `H` with two moons. |
| `solaris_system.json` | Solaris (*Solaris*) | A single ocean world, class `O`, around a binary pair — approximates the novel's living-ocean planet. |
| `trantor_system.json` | Trantor (*Foundation*) | Sun-like star, class `M` capital world, `age: "old"`, `intelligent_life: true`. |

None of these are meant to be scientifically exact reproductions — they map
each world onto the closest-fitting planet class this generator supports
(see `JSON.md`'s class table) and fill in the rest (star mass/radius within
its type, exact orbital distances, atmosphere specifics, names, flavor text)
normally. A couple of notable limitations worth knowing about:

- **Binary companions aren't independently typeable.** `star_type` only
  pins the primary; a binary system's secondary star always gets a randomly
  generated mass fraction of the primary. Solaris's actual red-star/blue-star
  pairing isn't achievable exactly for this reason.
- **Moon counts are best-effort.** A planet only has so much stable orbital
  room; `"moons": N` may generate fewer than `N` if there isn't room for all
  of them (see e.g. `krypton_system.json`, which asks for 2 but may get fewer).
