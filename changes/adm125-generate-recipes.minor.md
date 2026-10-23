**Headline:** A Generate form can be saved as a named recipe and opened again from the Start screen.

### Added
- **Generate recipes (ADM.125).** The Review step of every task that can be run again has a Save as recipe box. The Start screen lists your recipes: opening one shows the Review step filled with the saved settings (nothing runs until you press Run), and each can be renamed or deleted. Recipes belong to the admin who saved them; no other admin sees, opens, renames or deletes them. Wipe and Start a new galaxy cannot be saved. A name can be used once per admin, up to 80 characters, and an admin keeps at most 100.
- **Control schema v17 adds `generate_recipes` and `generate_recipe_fields`** (a recipe's admin, name and times, and its form one row per field). The control database has no Alembic chain, so the table arrives with `CREATE TABLE IF NOT EXISTS` on the next `planetgen.cli.migrate` run like the other control tables.
