**Headline:** A body dropped for the bar on habitable worlds no longer pulls tides on the bodies that stay.

### Changed

- GEN.1604: `generation/system._drop_barred_bodies` ran after the tide pass, so the tides, the heating and the effect-ledger rows of the bodies that stay still counted a planet or moon that was no longer in the system (on the 11-12 generator seed 433 of a bar-habitable system had 40 such rows). It now runs first, so the tide pass and the judgement of what the config asked for see the bodies that are left. Stored ocean tides, heating and the ledger change for a system that drops a body.

generation-output: changed
generation-affects: fill
