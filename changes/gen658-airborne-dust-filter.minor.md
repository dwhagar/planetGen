**Headline:** A world's drawn dust haze now decides whether a visitor should wear a dust filter.

### Changed
- **Airborne dust feeds the dust filter (GEN.658).** The stored dust optical depth (`tau_dust`) and the air's scale height give the dust in the air at the ground (`aerosols.dust_loading_mg_m3`, in mg/m3). From 0.1 mg/m3 (the NASA limit for Mars dust over 30 days) a dust filter is recommended, and from 10 mg/m3 it is required; a world with air and a dusty dry surface now lists it, and a sealed suit reads it as dust-tight seals. A body under 10 kPa still needs the filter or seals as before. The toxicity half waits for the surface composition to be stored. Derived from stored values, so there is no epoch or migration; `planetgen rescore-habitability` updates old galaxies.
