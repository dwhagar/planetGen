**Headline:** The info panel's Colour line now says which stored values gave a body's colour and whether it is a partial estimate or the full account.

### Changed

- `bodylook.colour_account` builds the line for a star (full when its colour is a model atmosphere at its temperature, partial from its spectral class alone), a solid world (ground colour from its mean surface temperature, sea colour from its ocean class, limb glow from its air under its star) and a giant (cloud-top colour from its temperature), each naming a default where an input is not stored and what is not used yet (surface composition, clouds and haze). `look_note` puts it ahead of the spin and lock notes. generation-output: unchanged
