**Headline:** The world classes and the life classes now judge oxygen by its pressure and not by its share of the air, so a thick air with a breathable oxygen pressure and a thin air with plenty of oxygen are no longer misjudged (GEN.1321).

### Changed

- The temperate world classes (M, H, O, P, R) take 16 to 35.5 kPa of oxygen, the marginal thin-air world (L) 3 to 16 kPa and the anoxic world (F) under 3 kPa, the Earth-air equivalents of the old 16, 35 and 3 percent; the fit tolerance (CF7) is a pressure tolerance for oxygen (stand-ins OP1, OP2 in section 6.196).
- The carbon dioxide edge of the oxygen-water classes is a pressure too: 1 kPa (class L from it, classes M, H, O, P and R below it); classes F and G keep the carbon dioxide share. `tests/class_fit_cli.py --compare REV` counts the bodies that change class code against an older fit table.
- Life class 2 needs under 3 kPa of free oxygen and class 4 at least 100 Pa (`o2_partial_pa`, a new quantity of the class and life-class rules), with the help text and the rule tables following.
- Oxygen as redox chemistry (the reducing-air test of the age needs) stays on the share.
