**Headline:** The class filters in search show each class's description again, every class has a marker colour on the system map, and the class words and docs match the rules.

### Fixed

- The planet and moon class filters in search read a description that was removed when class text moved to its own table (GEN.1212), so every class tooltip was empty. They show the class's own sentence now.
- The system map drew classes Q, R, U, W, X, Y and Z in the grey default colour. Each now has its own colour, and a test checks that every class letter has one and that no two share it.
- Class O's text said water covers "most" of the surface, though the rule is 95 percent or more (a global ocean with at most islands); it says "nearly all" now. Class X's text said "a thin air or none", though the rule allows any air above the ice lid (GEN.1737); it says "an air above it or none".
- The system file reference said Q is a reserved letter and cannot be used, that O starts at 80 percent water and M ends at 80 percent, and that class T takes a rocky body by its envelope share. It now matches the rules (Q is the hot world, O is 95 percent or more, M is 20 to 95 percent, T is hydrogen 30 percent or more of the air under 10 Earth masses).

generation-output: unchanged
generation-affects: none
