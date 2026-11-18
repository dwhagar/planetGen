**Headline:** A world that rolled the other-branch oddity stops looking for its second climate once two starts end on the answer it already has, so a dense warm planet no longer costs a dozen extra air solves.

### Changed

- W20 (`air_branches.other_branch`) solved the air again from `COLD_STARTS_K` (four starts) or `WARM_STARTS_K` (three) until one ended on a distinct answer. Of 132 worlds that rolled W20 in the redraw census systems 0 to 79, 129 had every start end on the first answer and one found a second branch (from the third start, after two that did not agree with the first answer). The search now stops when `AGREEING_STARTS_TO_STOP` (2) starts have ended less than `MIN_BRANCH_GAP_K` and `MIN_ICE_GAP` from the first answer, a start that ends elsewhere not counting. The stand-in row NEB8 names the new bound (tier T6, tweak level L1).
- Why now: the Train 45 head took 192 s for `world-requested` system 10 (7 s on main) because the request redraws a habitable world 68 times and half of those rolled W20, each paying three or four 2.5 to 3 s solves of a 2 bar air; 139 of the first 182 s were in `other_branch`. The test limit of 180 s is unchanged.
- Output: a world that rolled W20, whose first two starts agree with its first answer and whose later start would have found a distinct answer, now keeps its first answer. None of the 132 sampled does. Every other world, and every world that did not roll W20, is unchanged.

generation-output: changed
generation-affects: fill
