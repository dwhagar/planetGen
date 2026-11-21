**Headline:** A cold sulfuric acid cloud is no longer taken liquid down to 185 K: it stands only where its droplets are above their freezing point.

### Changed
- GEN.1644: the acid cloud step (`sulfur_cloud_set.acid_cloud_of_column`) now reads the freezing curve of the sulfuric acid hydrates (`acid_freezing`, [R4877: Figure 1]): a level whose droplets are at or below the freezing point of their strength is not part of the liquid cloud, so its base is the warmest liquid level and an air whose acid is liquid nowhere has no acid cloud. In the 5,000-system census one rocky body in 13,791 loses its cold acid cloud (a class E world with the ground at 198 K had a 23.5-deep one) and 20 more only have a colder top; the solid cloud and supercooled droplets are not modelled. Stand-in SV7 of section 6.228 is rewritten from a floor to this rule (T5, L0).

generation-output: changed
generation-affects: fill
