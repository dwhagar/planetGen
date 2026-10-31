# tests/test_autolink.py

"""
The auto-linker (`web/lib/autolink.py`): the first two mentions of a term on a page link to its help page, links already there count,
headings and existing links are left alone, table cells are linked, and the terms come from the parts of the program that registered them.
"""

from planetgen.web.lib import autolink, help_registry
from planetgen.web.lib.help_registry import Term

URL = {"4c": "/help/life/classes/4", "2": "/help/life/classes/2", "Class 3 family": "/help/life/classes/3"}
TERMS = [Term("4c", "x", {}), Term("2", "x", {}, "Class "), Term("Class 3 family", "x", {})]


def run(body, own=None):
    return autolink.link_page(f"<html><header>Class 4c</header><main>{body}</main></html>", TERMS, lambda t: URL[t.text], own)


def test_the_first_two_mentions_are_linked_and_the_rest_are_not():
    html = run("<p>Class 4c, 4c and 4c life.</p>")
    assert html.count('<a href="/help/life/classes/4">4c</a>') == 2
    assert html.endswith("and 4c life.</p></main></html>")


def test_a_link_already_there_counts_toward_the_two():
    one = run('<p><a href="/help/life/classes/4">life</a> and 4c, 4c.</p>')
    assert one.count('href="/help/life/classes/4"') == 2 and one.rstrip().count(">4c</a>") == 1
    two = run('<p><a href="/help/life/classes/4">a</a> <a href="/help/life/classes/4">b</a> 4c</p>')
    assert two.count("<a ") == 2 and "4c</p>" in two


def test_headings_links_and_the_site_header_are_left_alone_and_table_cells_are_linked():
    html = run('<h2>Class 4c</h2><p><a href="/other">Class 4c</a></p><table><tr><td>4c</td></tr></table>')
    assert "<h2>Class 4c</h2>" in html and '<a href="/other">Class 4c</a>' in html
    assert '<td><a href="/help/life/classes/4">4c</a></td>' in html
    assert html.startswith("<html><header>Class 4c</header>")


def test_a_bare_number_links_only_after_class_and_the_prefix_stays_outside_the_link():
    html = run("<p>Class 2 has 2 kinds in 2b.</p>")
    assert 'Class <a href="/help/life/classes/2">2</a> has 2 kinds in 2b.' in html


def test_whole_words_only_and_the_longer_term_wins():
    html = run("<p>x4c 4cc Class 3 family</p>")
    assert "x4c 4cc" in html and '<a href="/help/life/classes/3">Class 3 family</a>' in html


def test_a_page_does_not_link_to_itself_and_entities_survive():
    html = run("<p>4c &amp; 4c</p>", own="/help/life/classes/4")
    assert "<a " not in html and "4c &amp; 4c" in html


def test_a_page_without_main_or_terms_is_returned_as_it_came():
    assert autolink.link_page("<p>4c</p>", TERMS, lambda t: URL[t.text]) == "<p>4c</p>"
    assert autolink.link_page("<main>4c</main>", [], lambda t: "") == "<main>4c</main>"


def test_the_registered_life_terms_link_a_real_class_id():
    terms = help_registry.terms()
    html = autolink.link_page("<main><p>Class 4c life.</p></main>", terms, lambda t: "/p/" + str(t.values.get("number")))
    assert '<a href="/p/4">4c</a>' in html


def test_a_comment_pair_keeps_the_linker_out_of_the_text_between_them():
    html = run("<p>4c <!--autolink:off-->4c<!--autolink:on--> 4c 4c</p>")
    assert "<!--autolink:off-->4c<!--autolink:on-->" in html
    assert html.count('">4c</a>') == 2
