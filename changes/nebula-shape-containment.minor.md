### Changed
- A nebula holds a star system, phenomenon or smaller nebula only when the point lies inside its shape (GEN.75), not anywhere in its bounding sphere. Supernova remnants stay spheres. Containment is recomputed when a sector or nebula is generated; to refresh existing data, regenerate the affected sectors.
