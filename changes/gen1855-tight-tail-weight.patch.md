**Headline:** Algol binaries now reach the catalogue's lower limit of 0.1% of a disk of constant star formation.

### Fixed
- **GEN.1855: the Algol share sat under its lower limit.** The tight-tail weight (the factor on the density of periods under a day for stars lighter than the massive limit, `system_types.TUNING_KNOBS["tight_tail_weight"]`) is 8.0, up from 6.0. Over 1.5 million draws of a disk of constant star formation (ages 0.05 to 10 Gy) the share of systems with a semi-detached pair is 0.0881% at 6.0, 0.0952% at 7.0 and 0.1002% at 8.0; nothing had drifted since GEN.270 measured 0.09%. A tighter period tail puts more pairs under a day in every system of a lighter star, so a tree's pairs and every pair under a day shift: the contact share of FGK dwarfs rises from 0.165% to 0.181%.

generation-output: changed
generation-affects: scatter, fill
