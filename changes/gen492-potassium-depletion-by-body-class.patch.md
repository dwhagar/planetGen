**Headline:** A sourced potassium depletion by body class is built and tested, ready for the generator; no world changes yet.

### Changed

- **Potassium depletion (GEN.492, lock-free step).** `potassium_depletion` holds the K/U of the bulk silicate Earth, Mars, Mercury, the Moon, Vesta and the angrite parent body (Dauphas et al. 2022 Table 1) over CI chondrites, and draws a body's depletion from its class: planets of 0.03 Earth masses or more centre on 0.19, smaller bodies on 0.015. Nothing calls it yet; `radiogenic` still draws log-uniform 0.03 to 1 until the wiring step. Registry rows KD1 to KD4.
