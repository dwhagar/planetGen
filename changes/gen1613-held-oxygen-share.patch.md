**Headline:** A deep, hot air whose mantle asks for more oxygen than the air has pressure now settles instead of keeping an unsettled air.

### Fixed
- **The gas equilibrium never settled for a held oxygen fugacity over the total pressure (GEN.1613).** With the 2,000 K cap gone (GEN.1580) an oxidised mantle under a 4,000 K surface and 4,000 to 9,000 bar asked for 10^4.1 bar of oxygen in 9,255 bar, which an ideal gas cannot have, so the solve never settled (more iterations did not help; the cause was this, not the depth of the steam). `gas_equilibrium.solve` takes `lower_held_o2`; the rocky air (`mixture_rocky`) sets it, and a solve that failed at a held value over `MAX_HELD_O2_SHARE` (0.9) of the pressure is made again at that share, with one normal log line and `Equilibrium.clamped_o2_log_fugacity` holding the value asked for. A solve that settles is not touched. Stand-in HO1 (tier T3, tweak L3) in the registry; the real-gas fugacity coefficients that would let the mantle's value stand are their own item.

generation-output: changed
generation-affects: fill
