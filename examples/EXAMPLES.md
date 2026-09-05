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
| `van_maanens_star_system.json` | Van Maanen's Star (real) | A real, isolated white dwarf; covers the white-dwarf (`VII`) age/heliosphere path. |
| `zeta_ophiuchi_system.json` | Zeta Ophiuchi (real) | A real O-type main-sequence runaway star; covers the hot O/B dwarf wind path. |
| `vy_canis_majoris_system.json` | VY Canis Majoris (real) | One of the largest known stars; covers the extreme hypergiant (`0`) wind path. |
| `procyon_system.json` | Procyon (real) | A real F-type subgiant; covers the subgiant (`IV`) evolved-star age path. |

These last four exist mainly for **test coverage**: the original seven only
ever exercise Yerkes classes `V` (main sequence), `III` (giant), and `IB`
(supergiant) between them, missing exactly the white dwarf, hypergiant, and
O/B dwarf paths where real bugs were found and fixed (see git history for
`stellarObjects/physical_constants.py` and `starData.py`). `tests/test_examples.py`
runs every file in this directory automatically, so adding a new one here
adds it to the regression suite for free.

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
