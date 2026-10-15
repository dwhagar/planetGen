"""Shared helpers every planetgen package can use.

Logging and settings: `log` (the one console and file logging channel), `logpaths`
(where logs go; standard library only, so the installers can run it),
`settings` (the Pydantic model of every `config.json` option) and `stagemark` (the
"Stage 4 of 12" wording). Randomness: `draw` (every random draw generation
makes), `random` and `poisson`. Numbers and text: `format` (every measurement as text),
`ladders` (the unit ladders behind it), `checks` (`finite_domain`) and
`serialization` (the allowlist `to_dict` and `from_dict` helpers)."""
