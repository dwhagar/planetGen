**Headline:** The oxygen and carbon dioxide rules were walked end to end after the move to pressures: the class sentences, the class-fit report, the stand-in rows and the design notes now say kPa where the rules read kPa, and a test keeps it so (GEN.1331).

### Changed

- The sentences of classes L, F and M (planet_classes.json) give oxygen and carbon dioxide in kPa, as the rules read them.
- The class-fit report prints the oxygen and carbon dioxide pressure of a body beside the share.
- Stand-in rows CL3, CL4 and CF7 and the class note in condition-based-planet-classes.md say which edges are pressures and which stay shares.
- New `physics/data/partial_pressure_rules.json` lists every converted rule, the rules that stay shares with the reason (classes F and G, life class 6, the reducing-air and carbon-dioxide-air age rows) and the pressure constants; `tests/test_partial_pressure_rules.py` fails when a share-based oxygen or carbon dioxide threshold appears in a converted rule or in a rule the list does not name.

generation-output: unchanged
