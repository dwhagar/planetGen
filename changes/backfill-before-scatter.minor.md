**Headline:** The galaxy now fills in dimmer stars around the filled sectors before the big scatters, and the luminous-star default is 18,000 solar luminosities.

### Changed
- **Backfill runs before the scatters.** The galaxy run now goes sectors, link, clear, backfill, mass scatter, luminosity scatter, phenomena, then the rest. The scatter limits are recorded once, before the backfill, so the backfill is bounded by the same floor and mass limit the scatters then fill to, and the two hand-offs no longer disagree. The Generate page and the job runner list the stages in the same order. The order was picked from tiny-galaxy runs of every order of the post-planning steps (reports/generation-order-2026-11-13/report.md).
- **Luminous-star default doubled.** The default luminosity floor is 18,000 solar luminosities (was 9,000). In the default galaxy that is about 5,983 luminosity-pass stars instead of 24,794.
- generation-affects: scatter, fill
- generation-output: changed
