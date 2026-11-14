**Headline:** A request for a habitable world now asks for land and sea and for air of 50 to 250 kPa.

### Changed
- **A habitable-world request asks for 20 to 95 percent water cover and 50 to 250 kPa of air (GEN.1702).** A planet or moon meets the request (`habitable_request.meets`) only when its surface is temperate, its water is a surface ocean, the liquid ocean covers 20 to 95 percent of the surface (`tuning.REQUEST_WATER_COVER`, so there is land and there is sea) and its total air pressure is 50 to 250 kPa (`tuning.REQUEST_PRESSURE_PA`, replacing the 10 kPa floor `REQUEST_MIN_PRESSURE_PA`). Boss's card "Yes, add land and pressure", 11:18 PM EST 2026-11-13. A request still only raises the odds and never forces (GEN.1130); a system that bars habitable worlds keeps out a body that meets the new conditions or is rated habitable. Stand-in rows CC11 (cover) and CC12 (pressure) of the registry are rewritten. Evidence: 109 Sun-like systems with the request gave 0 Class M, and the kept planets were 100 percent ocean (reports/class-m-blue-levers-2026-11-14.md, section 3).
- **Reuse:** the existing `habitable_request.meets` test and `tuning.REQUEST_*` constants; no new code path and no shim.

generation-output: changed
generation-affects: fill
