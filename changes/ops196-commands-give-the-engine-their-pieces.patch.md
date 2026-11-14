**Headline:** A guard test now keeps every queue command on the engine's loop, so no new command can go silent by writing a loop of its own.

### Changed

- OPS.196 (first part): `test_commands_have_pieces` fails for a registered command that gives the engine neither `pieces` and `do` nor a place in its list of commands waiting for the family item that converts them (OPS.198 to OPS.206). The list only gets shorter: a converted command still listed fails, and so does a name that is not a command. The second part (each command run on the small galaxy with its longest quiet stretch held to a budget, and the `whole` flag deleted) comes with the last family.

generation-output: unchanged
generation-affects: none
