**Headline:** a pure loss term for the heavy gases of an air (the replacement for the stepped pressure caps) is built and tested, and nothing calls it yet.

### Changed

- New `physics/heavy_gas_loss.py` (GEN.1558 prebuild, stacked on nothing, wired nowhere): `lifted_mass_kg` sums `thermosphere.energy_limited_mass_loss_kg_s` over the same log-spaced XUV history and decay as `thermosphere.kept_shares`; `kept_share` turns it into the share of an air that stays, `exp(-efficiency * open_share * lifted / air_mass)`, with the open share of `ion_escape`. It reuses both; `atmosphere_retention.energy_limited_loss_fraction` is not used.
- New test `test_heavy_gas_loss.py` (6 tests): bounds, monotonic in flux, mass, efficiency and field, and one efficiency (0.007) that lets an Earth with a field at ten radii keep its air while an open Mars and an open Moon-sized body lose theirs.
- The efficiency is an argument, not a default: the build must register it as a stand-in (GEN.725) and calibrate it. The measured window for one value is 0.005 to 0.015.

generation-output: unchanged
generation-affects: none
