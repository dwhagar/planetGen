**Headline:** A world's class is now read once, from the finished body after its life pass, as Boss ruled (no class until the evolution is finished); on 60 generated systems no class changes.

### Changed

- GEN.1745 slice 3 (class last): `life._rebuild` and `life._remake_ground_if_thawed` no longer read the class in the middle of the life pass, and `StarSystem._finish_planet` and `_finish_moon` prepare the class (the body type and the roundness ratio, `planets.prepare_class`) without classifying; `life.apply_life_data` reads the class once, at its end, from the body as the life pass left it (the air, the water and the score after the life, the trace gases and the stage it reached). A small body that has air or life from its life pass is now classed by that finished body, never by the one it had before (GEN.1752). Nothing in the generator reads a body's class between those points. Measured on 60 generated systems (1,247 bodies) against main: no body changes class, variant or match.

generation-output: changed
generation-affects: fill
