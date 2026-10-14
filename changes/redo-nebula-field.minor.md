### Added
- `planetgen plan --redo-nebula-field [dry-run|apply]`: brings an existing galaxy's stored molecular clouds in line with the current keep rates. The dry run (the default) counts the stored field clouds the rule now rejects, per class, and changes nothing; `apply` deletes them, refreshes containment in the sectors they reached and draws the `nebula_field` table again (GEN.205).
