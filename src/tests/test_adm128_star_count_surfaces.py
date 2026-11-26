# tests/test_adm128_star_count_surfaces.py

"""
ADM.128: the one-off system page and `POST /api/systems` take `stars` and `max_stars` and refuse what `planetgen system --stars
--max-stars` refuses, in the command line's words (`run_system.stars_problem`, `max_stars_problem`): a count under 1, over the cap,
over 64, a cap outside 1 to 64, one star with a binary feature, and two or more stars with `binary_system` off.
"""

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.web import generate_page, system_page
from tests.test_api import admin_client, client, default_admin_client, first_admin_password  # noqa: F401 -- fixtures

# (the API body, the page's form, the command line's options, what all three say)
REFUSALS = {
    "stars under 1": ({"stars": 0}, {"stars": "0"}, ["--stars", "0"], "stars 0 must be a whole number of at least 1."),
    "stars over the cap": ({"stars": 13}, {"stars": "13"}, ["--stars", "13"], "more than the cap of 12 stars"),
    "stars over a lowered cap": ({"stars": 9, "max_stars": 8}, {"stars": "9", "max_stars": "8"}, ["--stars", "9", "--max-stars", "8"],
                                 "more than the cap of 8 stars"),
    "stars over 64": ({"stars": 65, "max_stars": 64}, {"stars": "65", "max_stars": "64"}, ["--stars", "65", "--max-stars", "64"],
                      "more than any system may hold (64)"),
    "cap under 1": ({"max_stars": 0}, {"max_stars": "0"}, ["--max-stars", "0"], "max_stars 0 must be a whole number from 1 to 64."),
    "cap over 64": ({"max_stars": 65}, {"max_stars": "65"}, ["--max-stars", "65"], "max_stars 65 must be a whole number from 1 to 64."),
    "one star, binary": ({"stars": 1, "binary_system": True}, {"stars": "1", "binary_system": "yes"}, ["--stars", "1", "+binary_system"],
                         "single star, which +binary_system contradicts"),
    "one star, wide binary": ({"stars": 1, "wide_binary": True}, {"stars": "1", "wide_binary": "yes"}, ["--stars", "1", "+wide_binary"],
                              "no companion for +wide_binary"),
    "two stars, binary off": ({"stars": 2, "binary_system": False}, {"stars": "2", "binary_system": "no"}, ["--stars", "2", "-binary_system"],
                              "more than one star, which -binary_system contradicts"),
    "seven stars, binary off": ({"stars": 7, "binary_system": False}, {"stars": "7", "binary_system": "no"}, ["--stars", "7", "-binary_system"],
                                "more than one star, which -binary_system contradicts"),
}


@pytest.mark.parametrize("case", sorted(REFUSALS))
def test_the_api_refuses_what_the_command_line_refuses_in_its_words(case, admin_client):
    body, _form, _argv, words = REFUSALS[case]
    response = admin_client.post("/api/systems", json=body)
    assert response.status_code == 400, response.get_json()
    assert words in response.get_json()["error"]


@pytest.mark.parametrize("case", sorted(REFUSALS))
def test_a_regeneration_recipe_refuses_the_same_cases(case):
    from planetgen.api.schemas import RegenerateRecipe
    from planetgen.api.common import ApiError
    from planetgen.api.schemas import parse_body

    body, _form, _argv, words = REFUSALS[case]
    with pytest.raises(ApiError, match=words.replace("(", r"\(").replace(")", r"\)").replace("+", r"\+")):
        parse_body(RegenerateRecipe, body)


@pytest.mark.parametrize("case", sorted(REFUSALS))
def test_the_page_refuses_what_the_command_line_refuses_in_its_words(case):
    _body, form, _argv, words = REFUSALS[case]
    with pytest.raises(generate_page.FormError) as caught:
        system_page.system_request(form)
    assert words in str(caught.value) or words[0].upper() + words[1:] in str(caught.value)


@pytest.mark.parametrize("case", sorted(REFUSALS))
def test_the_command_line_says_the_same(case, capsys):
    _body, _form, argv, words = REFUSALS[case]
    parser, parsers = generate_cli.build_parser()
    args = parser.parse_args(["system"] + argv)
    with pytest.raises(SystemExit):
        generate_cli.validate_system_args(args, parsers["system"])
    assert words in capsys.readouterr().err


@pytest.mark.parametrize("body, form", [
    ({"stars": 3, "max_stars": 5}, {"stars": "3", "max_stars": "5"}),
    ({"stars": 1, "binary_system": False}, {"stars": "1", "binary_system": "no"}),
    ({"stars": 3, "binary_system": True}, {"stars": "3", "binary_system": "yes"}),
    ({"stars": 64, "max_stars": 64}, {"stars": "64", "max_stars": "64"}),
    ({"max_stars": 1}, {"max_stars": "1"}),
])
def test_a_count_that_fits_is_accepted_by_both(body, form):
    from planetgen.api.schemas import SystemCreate, parse_body

    assert parse_body(SystemCreate, body).model_dump(exclude_unset=True) == body
    spec = system_page.system_request(form)
    for key in ("stars", "max_stars"):
        if key in form:
            assert f"--{key.replace('_', '-')}={form[key]}" in spec["argv"]


def test_a_count_must_be_a_whole_number_in_the_api():
    from planetgen.api.common import ApiError
    from planetgen.api.schemas import SystemCreate, parse_body

    for body in ({"stars": "3"}, {"stars": 2.5}, {"stars": True}, {"max_stars": "9"}):
        with pytest.raises(ApiError):
            parse_body(SystemCreate, body)


def test_the_api_builds_a_system_with_the_forced_count(admin_client):
    response = admin_client.post("/api/systems", json={"stars": 3, "max_stars": 5, "planets": False})
    assert response.status_code == 201, response.get_json()
