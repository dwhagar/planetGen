**Headline:** The known-issue ship rule now says a late output-moving add-on to a train re-runs the pins and the golden seed on the final head.

### Changed
- **Late add-on rule** (TEST.725): `docs/plan/lanes.md` and the Known issues note of the train PR body (`scripts/train.py pr-body`) say that an output-moving fix added after the gate passed (Train 42's GEN.1774 moved the civilized pins after the gate) is not shipped on the gate's word; the pins and the golden seed run again on the final head. A fix found in the gate and put into the train before its final gate needs no extra run.
