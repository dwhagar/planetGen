# Retention and lava worlds: the check of rows d1, i6 and j1 (GEN.511)

Status: built 2026-10-23 by lane M (builder B3). Rows d1 (the cosmic shoreline), i6 (a molten dayside keeps volatiles) and j1 (the retained fraction after a runaway) of the atmosphere recipes (`scripts/atmosphere_recipes_data.py`); tests in `src/tests/test_tuned_draws.py`.

## 1. d1: the shoreline

Row d1 already draws the published logistic: the chance a planet keeps an atmosphere is `1 / (1 + 10^-(z))` with `z = (log10 f0 - log10 S) / 0.247`, where `log10 f0 = 2.62` at the Earth's escape velocity and a solar luminosity [R0364], the escape velocity enters with `p = 5.9` and the luminosity with `q = 1.17`. The builder found nothing to retune; it pinned the row with a test of the four cases (half at the line, one width above, the escape velocity slope, the luminosity slope) so a later edit cannot move it silently.

## 2. i6: the molten-dayside share

The row's evidence is 2 of 2 molten-dayside ultra-short-period planets that keep volatiles (N = 2). The bracket was 0.1 to 0.9 with a note that said 0.29 to 0.99. The note was the correct statistic: the 95 percent interval of Beta(3, 1) (a flat prior and 2 successes in 2) is 0.29 to 0.99. The bracket now matches it. With N = 2 the row stays partly verified, and no target exists (a target needs a count with a spread to compare against).

## 3. j1: the retained fraction

The runaway gate is [R2446] and is not tuned. The retained fraction was drawn log-uniform from 0.1 to 1. Of the two bare M-dwarf planets of the sample (TRAPPIST-1 b and c), none kept an atmosphere, so a lower bound of 0.1 forbids what is observed. The draw is now uniform on 0 to 1, which allows a bare planet and every retained fraction up to all of it.

## 4. What it does not do

- The recipes are not yet read by the generator, so no generated value changes and no epoch is taken.
- Whether TRAPPIST-1 b and c are bare is itself uncertain (JWST thermal measurements allow a thin atmosphere), so the row is not tuned to a share; the bracket only keeps the observations possible.
