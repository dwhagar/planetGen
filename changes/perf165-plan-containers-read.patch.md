**Headline:** A galaxy plan no longer reads every stored cluster back just to build the containers of its globular ones.

### Changed
- The plan's container step reads only the King-profile clusters (157 of 37,402 at the default size) instead of all of them, which takes about 1 to 1.5 s off a plan of a default galaxy. The containers are the same.
- `docs/design/generation-performance.md` records why the component insert, the second `plan_shape` and the start-up math check were left as they are.
