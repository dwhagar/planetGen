**Headline:** A hot air keeps its equilibrium nitrogen dioxide and nitrous oxide when the trace gases are drawn.

### Fixed
- **`trace_air` removed the equilibrium NO2 and N2O of a hot air (GEN.1614).** Before drawing the trace gases it took every trace gas out of the air, so the NO2 and N2O that the equilibrium holds in a hot air (4,000 K and more, once GEN.1580 lets a surface get there) went with the drawn ones. It now takes off each trace gas only what a draw can have given it, the top of the gas's design bracket (`trace_air.TRACE_CEILING`: N2O 1e-6, PH3 2e-8, NO2 1e-8, CH3Cl 1e-9 of the air), and keeps the rest.

generation-output: changed
generation-affects: fill
