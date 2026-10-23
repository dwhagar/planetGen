# planetgen/admin/recipes.py

"""
Generate recipes (ADM.125): a named, saved Generate form, kept per admin in the control database
(`generate_recipes` and its `generate_recipe_fields`, `control_schema.sql`).

A recipe is what "Run again" already keeps for a job, a draft: the fields of the task's form as `[name, value]`
pairs, the first of which names the action, kept one `generate_recipe_fields` row per pair. Saving one is done from the Review step; the Start screen lists them,
a recipe opens the task's Review step filled with its draft (nothing runs until Run is pressed there), and it can be
renamed or deleted. Each is the admin's own: every function takes the admin's account name and only ever reads or
changes that account's rows, so one admin cannot see, run, rename or delete another's.

This module imports no Flask; the web layer (`web/generate_flow.py`) calls it.
"""

import re

from planetgen.db import store

NAME_MAX = 80
"""int: The longest recipe name, in characters (the column's width)."""

MAX_RECIPES = 100
"""int: How many recipes one admin may keep."""

MAX_FIELDS = 300
"""int: The most fields a draft may hold (the largest Generate form has about a hundred)."""

FIELD_NAME_MAX = 100
"""int: The longest field name a draft may hold, in characters (the column's width)."""

MAX_VALUE = 4000
"""int: The longest field value a draft may hold, in characters (the column's width)."""


class RecipeError(ValueError):
    """A recipe that cannot be saved or changed; the message says why, in words for the admin."""


def clean_name(name):
    """`name` with its white space tidied (runs of it become one space).

    Raises:
        RecipeError: The name is empty or longer than `NAME_MAX`.
    """
    tidy = re.sub(r"\s+", " ", str(name or "")).strip()
    if not tidy:
        raise RecipeError("Give the recipe a name.")
    if len(tidy) > NAME_MAX:
        raise RecipeError(f"A recipe's name is at most {NAME_MAX} characters.")
    return tidy


def _clean_draft(draft):
    """`draft` as a list of `[name, value]` string pairs.

    Raises:
        RecipeError: It holds no fields, too many, or one too long.
    """
    pairs = [[str(name), str(value)] for name, value in draft]
    if not pairs:
        raise RecipeError("This form cannot be saved as a recipe.")
    if len(pairs) > MAX_FIELDS or any(len(a) > FIELD_NAME_MAX or len(b) > MAX_VALUE for a, b in pairs):
        raise RecipeError("This form is too large to save as a recipe.")
    return pairs


def _admin_id(conn, username):
    row = conn.execute("SELECT id FROM admin_users WHERE username = ?", (username,)).fetchone()
    if row is None:
        raise RecipeError("Only a signed-in admin can keep recipes.")
    return row["id"]


def _recipes(conn, where, params):
    """The recipes `where` (a condition on the `generate_recipes` row `r` and its admin `a`) selects, by name, each with
    its draft read from its field rows."""
    rows = conn.execute(
        "SELECT r.id, r.name, r.created_at, r.updated_at FROM generate_recipes r "
        f"JOIN admin_users a ON a.id = r.admin_id WHERE {where} ORDER BY r.name, r.id", params).fetchall()
    found = []
    for row in rows:
        fields = conn.execute("SELECT field_name, field_value FROM generate_recipe_fields WHERE recipe_id = ? "
                              "ORDER BY position", (row["id"],)).fetchall()
        found.append({"id": row["id"], "name": row["name"], "created_at": row["created_at"],
                      "updated_at": row["updated_at"],
                      "draft": [[f["field_name"], f["field_value"]] for f in fields]})
    return found


def list_recipes(username, config=None):
    """The admin's recipes, by name: `[{id, name, draft, created_at, updated_at}]`."""
    conn = store.get_control_connection(config)
    try:
        return _recipes(conn, "a.username = ?", (username,))
    finally:
        conn.close()


def get_recipe(username, recipe_id, config=None):
    """The admin's recipe `recipe_id`, or `None` when there is no such recipe of theirs."""
    conn = store.get_control_connection(config)
    try:
        found = _recipes(conn, "a.username = ? AND r.id = ?", (username, recipe_id))
        return found[0] if found else None
    finally:
        conn.close()


def save_recipe(username, name, draft, config=None):
    """Keeps `draft` (`[name, value]` pairs) as the admin's recipe `name`; returns the recipe.

    Raises:
        RecipeError: The name or the draft is not acceptable, the admin has `MAX_RECIPES` already, or a recipe of
            theirs has the name (it is never replaced, so rename or delete the old one first).
    """
    name = clean_name(name)
    pairs = _clean_draft(draft)
    conn = store.get_control_connection(config)
    try:
        admin_id = _admin_id(conn, username)
        if conn.execute("SELECT COUNT(*) AS n FROM generate_recipes WHERE admin_id = ?",
                        (admin_id,)).fetchone()["n"] >= MAX_RECIPES:
            raise RecipeError(f"You have {MAX_RECIPES} recipes, the most there can be. Delete one first.")
        if conn.execute("SELECT 1 FROM generate_recipes WHERE admin_id = ? AND name = ?",
                        (admin_id, name)).fetchone():
            raise RecipeError(f"You already have a recipe named “{name}”. Rename or delete it first, or "
                              "choose another name.")
        recipe_id = conn.execute("INSERT INTO generate_recipes (admin_id, name) VALUES (?, ?)",
                                 (admin_id, name)).lastrowid
        for position, (field, value) in enumerate(pairs):
            conn.execute("INSERT INTO generate_recipe_fields (recipe_id, position, field_name, field_value) "
                         "VALUES (?, ?, ?, ?)", (recipe_id, position, field, value))
        conn.commit()
    finally:
        conn.close()
    return get_recipe(username, recipe_id, config)


def rename_recipe(username, recipe_id, name, config=None):
    """Renames the admin's recipe `recipe_id`; returns the recipe, or `None` when there is no such recipe of theirs.

    Raises:
        RecipeError: The name is not acceptable or another recipe of theirs has it.
    """
    name = clean_name(name)
    conn = store.get_control_connection(config)
    try:
        admin_id = _admin_id(conn, username)
        taken = conn.execute("SELECT id FROM generate_recipes WHERE admin_id = ? AND name = ? AND id <> ?",
                             (admin_id, name, recipe_id)).fetchone()
        if taken:
            raise RecipeError(f"You already have a recipe named “{name}”.")
        cur = conn.execute("UPDATE generate_recipes SET name = ? WHERE admin_id = ? AND id = ?",
                           (name, admin_id, recipe_id))
        conn.commit()
        if not cur.rowcount and conn.execute("SELECT 1 FROM generate_recipes WHERE admin_id = ? AND id = ?",
                                             (admin_id, recipe_id)).fetchone() is None:
            return None
    finally:
        conn.close()
    return get_recipe(username, recipe_id, config)


def delete_recipe(username, recipe_id, config=None):
    """Deletes the admin's recipe `recipe_id`; the deleted recipe's name, or `None` when there is no such recipe of
    theirs."""
    found = get_recipe(username, recipe_id, config)
    if found is None:
        return None
    conn = store.get_control_connection(config)
    try:
        conn.execute("DELETE r FROM generate_recipes r JOIN admin_users a ON a.id = r.admin_id "
                     "WHERE a.username = ? AND r.id = ?", (username, recipe_id))
        conn.commit()
    finally:
        conn.close()
    return found["name"]
