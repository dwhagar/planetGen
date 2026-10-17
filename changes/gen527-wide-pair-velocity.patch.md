**Headline:** The two stars of a wide pair now carry the same velocity when a sector is generated and when it is loaded back from the database.

### Changed
- **Every star of a system takes the system's velocity (GEN.527).** `place_system` used to give each star the rotation curve's tangent at its own place, plus what the system has beyond the curve. A wide pair's two stars sit far enough apart for the tangent to differ by 0.01 to 0.08 m/s, but the database stores one velocity per system, so a loaded pair came back with both stars on the system's vector and differed from the generated pair by that much (the cause of the occasional failure of `test_a_saved_runaway_system_loads_back_with_the_velocity_it_had`). Both paths now give each star the one system vector.
- The saved-system test is back to a relative bound of 1e-6, and a new test places wide pairs, clears the stars' own speeds as a load does, places them again from the stored vector, and requires identical star velocities.

TODO items resolved: GEN.527. Added: none.
