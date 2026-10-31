**Headline:** The property sheet now says why a giant planet shows no starlight row, and the fill audit places its sample systems in the galaxy so a star's lethal-event rate is checked.

### Changed
- **GEN.784.** `energy_flux_w_m2` and `flux_variation` are stored for every rocky or icy planet and every moon (12 sample systems: no exception), and left none for a giant by design (`insolation` and `habitability_world` score a giant's moons in its light but never the giant). The sheet now carries that rule with its reason and the open gap is removed.
- **GEN.785.** `lethal_event_rate_per_gyr` is set for every star from its place in the galaxy (`activity.set_galactic_hazard`) and stored with it; it was empty in the audit only because the sample systems were built with no galactic distance. The audit now places each sample 26,000 ly from the centre, as the galaxy does for every system it builds, and the open gap is removed. A test checks a placed star stores a positive rate.

Nothing a generator draws or stores changes.
generation-output: unchanged
