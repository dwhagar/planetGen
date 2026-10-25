**Headline:** The sulfate and sulfide balance that decides how much of the sulfur a volcano gives stays in the air is built as a module; nothing calls it yet.

### Added
- **Sulfur sink (GEN.815).** `physics/sulfur_balance.py` splits the degassed sulfur between the air, a calcium sulfate sink and an iron sulfide sink by the redox of the gas, limited by the crust the volcanism built, and leaves in the air only what the sinks have not yet taken (rain-out over liquid water, sulfation on a dry surface) and what they cannot hold. It replaces the stand-in that kept 1e-5 to 1e-2 of the sulfur in the air of a wet world once the surface-air step uses it. Five stand-ins are in the registry (AM45 to AM49). The module is not called by the generator yet, so no generated value changes.
