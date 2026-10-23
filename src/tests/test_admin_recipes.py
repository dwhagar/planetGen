# tests/test_admin_recipes.py

"""
Generate recipes (ADM.125, `planetgen/admin/recipes.py`): a named, saved Generate form kept per admin in the control
database. Needs the MySQL test server (the `mysql_config` fixture); skipped without it.
"""

import pytest

from planetgen.admin import recipes
from planetgen.db import store

DRAFT = [["action", "galaxy"], ["what", "grid"], ["mode", "shell"], ["shell_index", "40"]]


@pytest.fixture
def control_config(mysql_config, monkeypatch):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    for name in ("boss", "other"):
        conn.execute("INSERT INTO admin_users (username, password_hash) VALUES (?, 'x')", (name,))
    conn.commit()
    conn.close()
    monkeypatch.setenv(store.CONTROL_DB_ENV_VAR, mysql_config.database)
    return mysql_config


def test_a_saved_recipe_is_listed_with_its_draft(control_config):
    saved = recipes.save_recipe("boss", "  Shell   40 ", DRAFT, control_config)
    assert saved["name"] == "Shell 40" and saved["draft"] == DRAFT
    assert [r["name"] for r in recipes.list_recipes("boss", control_config)] == ["Shell 40"]
    assert recipes.get_recipe("boss", saved["id"], control_config)["draft"] == DRAFT


def test_the_longest_field_and_a_repeated_field_are_kept_in_order(control_config):
    draft = [["action", "galaxy"], ["n" * recipes.FIELD_NAME_MAX, "v" * recipes.MAX_VALUE], ["tag", "b"], ["tag", "a"]]
    assert recipes.save_recipe("boss", "long", draft, control_config)["draft"] == draft


def test_recipes_are_listed_by_name(control_config):
    for name in ("b", "c", "a"):
        recipes.save_recipe("boss", name, DRAFT, control_config)
    assert [r["name"] for r in recipes.list_recipes("boss", control_config)] == ["a", "b", "c"]


@pytest.mark.parametrize("name, message", [("", "name"), ("   ", "name"), (None, "name"),
                                           ("x" * (recipes.NAME_MAX + 1), "at most")])
def test_a_name_that_cannot_be_used_is_refused(control_config, name, message):
    with pytest.raises(recipes.RecipeError, match=message):
        recipes.save_recipe("boss", name, DRAFT, control_config)
    assert recipes.list_recipes("boss", control_config) == []


def test_the_longest_name_is_kept(control_config):
    name = "n" * recipes.NAME_MAX
    assert recipes.save_recipe("boss", name, DRAFT, control_config)["name"] == name


@pytest.mark.parametrize("draft", [[], [["a", "x" * (recipes.MAX_VALUE + 1)]], [["n" * (recipes.FIELD_NAME_MAX + 1), "v"]],
                                   [["a", "b"]] * (recipes.MAX_FIELDS + 1)])
def test_a_draft_that_is_empty_or_too_large_is_refused(control_config, draft):
    with pytest.raises(recipes.RecipeError):
        recipes.save_recipe("boss", "big", draft, control_config)


def test_a_name_in_use_is_never_replaced(control_config):
    first = recipes.save_recipe("boss", "same", DRAFT, control_config)
    with pytest.raises(recipes.RecipeError, match="already have"):
        recipes.save_recipe("boss", "same", [["action", "plan"]], control_config)
    assert recipes.get_recipe("boss", first["id"], control_config)["draft"] == DRAFT
    # another admin may use the name
    assert recipes.save_recipe("other", "same", DRAFT, control_config)["name"] == "same"


def test_an_admin_keeps_at_most_the_limit(control_config, monkeypatch):
    monkeypatch.setattr(recipes, "MAX_RECIPES", 2)
    recipes.save_recipe("boss", "one", DRAFT, control_config)
    recipes.save_recipe("boss", "two", DRAFT, control_config)
    with pytest.raises(recipes.RecipeError, match="Delete one"):
        recipes.save_recipe("boss", "three", DRAFT, control_config)
    recipes.save_recipe("other", "three", DRAFT, control_config)


def test_an_unknown_account_cannot_keep_recipes(control_config):
    with pytest.raises(recipes.RecipeError, match="signed-in admin"):
        recipes.save_recipe("nobody", "x", DRAFT, control_config)


def test_rename(control_config):
    one = recipes.save_recipe("boss", "one", DRAFT, control_config)
    recipes.save_recipe("boss", "two", DRAFT, control_config)
    assert recipes.rename_recipe("boss", one["id"], " first ", control_config)["name"] == "first"
    assert recipes.rename_recipe("boss", one["id"], "first", control_config)["name"] == "first"  # unchanged is fine
    with pytest.raises(recipes.RecipeError, match="already have"):
        recipes.rename_recipe("boss", one["id"], "two", control_config)
    with pytest.raises(recipes.RecipeError):
        recipes.rename_recipe("boss", one["id"], "", control_config)
    assert recipes.rename_recipe("boss", one["id"] + 999, "x", control_config) is None


def test_one_admin_cannot_see_rename_or_delete_anothers(control_config):
    mine = recipes.save_recipe("boss", "mine", DRAFT, control_config)
    assert recipes.list_recipes("other", control_config) == []
    assert recipes.get_recipe("other", mine["id"], control_config) is None
    assert recipes.rename_recipe("other", mine["id"], "stolen", control_config) is None
    assert recipes.delete_recipe("other", mine["id"], control_config) is None
    assert recipes.get_recipe("boss", mine["id"], control_config)["name"] == "mine"


def test_delete_returns_the_name_and_removes_the_recipe(control_config):
    mine = recipes.save_recipe("boss", "mine", DRAFT, control_config)
    assert recipes.delete_recipe("boss", mine["id"], control_config) == "mine"
    assert recipes.list_recipes("boss", control_config) == []
    assert recipes.delete_recipe("boss", mine["id"], control_config) is None


def test_deleting_an_admin_deletes_their_recipes(control_config):
    recipes.save_recipe("boss", "mine", DRAFT, control_config)
    conn = store.get_control_connection(control_config)
    try:
        conn.execute("DELETE FROM admin_users WHERE username = 'boss'")
        conn.commit()
        assert conn.execute("SELECT COUNT(*) AS n FROM generate_recipes").fetchone()["n"] == 0
        assert conn.execute("SELECT COUNT(*) AS n FROM generate_recipe_fields").fetchone()["n"] == 0
    finally:
        conn.close()
