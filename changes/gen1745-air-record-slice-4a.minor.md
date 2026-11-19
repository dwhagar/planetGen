**Headline:** A rocky world's air history now has a record: the steps with their temperature, pressure and gases, the mass extinctions and the life ages, in a short text form (nothing stores it yet)

### Added

- GEN.1745 slice 4a: `physics/air_record.py` assembles, for a rocky body with an air of 500 Pa or more, the recorded steps (at most 60: the history's schedule, the oxygen's rise, the producers' ages and each mass extinction) with each step's surface temperature, total pressure, partial pressure of every gas of 0.1 Pa or more and oxygen law, and the life ages; `pack` and `unpack` turn it to and from the JSON text a body's column will hold (slice 4b). The inputs of the final air solve are kept on the body's build record (`air_input`, not stored).

generation-output: unchanged
