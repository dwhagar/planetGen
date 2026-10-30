**Headline:** A world whose star is younger than the first step of its life schedule is no longer given a life class, so the map and the object pages no longer show life that has not started.

### Fixed

- GEN.1001: when a body's schedule puts abiogenesis after its star's age (cold, slow worlds round a young star), the life pass stores no class, schedule, stage or complexity for it, as for a world no class supports; before, it kept the class with no stage, and the class markers, the Has Class N life line and the class filters showed life that had not begun. This changes what is drawn for those worlds, so it takes a generator epoch (affects fill). The standing test that a body with a class has a stage no longer expects to fail.
