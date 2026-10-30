**Headline:** The generator can now tell where a world holds liquid sulfuric acid, the medium of acid-world life.

### Changed

- GEN.888: a new reader works out from a world's air whether it has a sulfuric acid cloud layer (like Venus's temperate clouds) or an acid pool on warm dry ground, and how strong the acid is. It needs sulfur dioxide, an oxidiser and a trace of water, keeps to the 283 to 373 K liquid window and 0.3 to 3 bar for clouds, and counts acid of 81 percent by weight or stronger. Nothing is stored and nothing changes on a generated world yet; the life class rules read it when the wiring lands with the next generation epoch.
- The air-profile helper accepts the temperature the air turns isothermal at, so a hot thick air gets a cold top from the absorbed starlight instead of a fixed fraction of its surface temperature.
