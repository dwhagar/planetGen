**Headline:** A system file's planet slot can ask for conditions, and a slot that cannot be read now fails in plain words.

### Changed
- **A planet slot names conditions (GEN.688).** The system file (command line and the Generate page's System file box) takes `"conditions"` in a slot: `kind`, `air` and a `[low, high]` range for a quantity such as `temperature_k` or `pressure_pa`. `"planet_class"` stays as a shorthand for the conditions of that class. The Generate page's example slot now shows conditions.
- **A bad slot is refused before anything is generated.** `run_system.slots_problem` says what is wrong (an unknown key, a type that is not planet or asteroid_belt, both a letter and conditions, a quantity no class reads, a range that is not two numbers, a negative moon count). The command line gives a usage error and the page shows it under the form.
