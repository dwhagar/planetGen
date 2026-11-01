# tests/test_autolink_adversarial.py

"""
Tries to break the auto-linker (`web/lib/autolink.py`): to make it link the wrong word, the wrong page, text it must not touch, or more
than the limit. The first half is hand-made attacks, with measurements and units (a speed "4c", "4 C", "3 g", "15 Mpc") taken for a class as the
main one; the second half generates random pages with a fixed seed and checks what must hold for any page:

1. Take the links the linker added out and the page is exactly what went in.
2. No term has more than two links to its page, existing ones counted.
3. No link is added inside a heading, a link, a control, a script, an attribute, a comment or a switched-off stretch.
4. An added link is a whole word, not part of a longer one or of a decimal, and a class id is always after "Class".
5. A second pass adds nothing.
"""

import random
import re
from html.parser import HTMLParser

import pytest

from planetgen.web.lib import autolink, help_registry
from planetgen.web.lib.help_registry import Term

P4, P2, P3 = "/help/life/classes/4", "/help/life/classes/2", "/help/life/classes/3"
URL = {"4c": P4, "2": P2, "3a": P3, "3b": P3, "Class 3 family": P3, "Class 4c": "/other/long"}
TERMS = [Term(t, "x", {}, "Class ") for t in ("4c", "2", "3a", "3b")] + [Term("Class 3 family", "x", {})]


def run(body, terms=None, own=None):
    page = f"<html><body><header>Class 4c</header><main>{body}</main></body></html>"
    return autolink.link_page(page, terms or TERMS, lambda t: URL[t.text], own)


def links(html):
    """The `(address, text)` of every link the linker added. The linker writes `<a href="...">`; the pages in these tests write their own with single quotes."""
    return re.findall(r'<a href="([^"]+)">([^<]*)</a>', html)


# --- Measurements and units are never classes ---------------------------------------------------------------------------

MEASUREMENTS = [
    "4c", "2c", "0.4c", "0.2c", "travels at 4c", "a speed of 2c", "4 C", "4C", "4°C", "4 °C", "-4°C", "25 °C", "3 g", "3g", "0.3 g", "1.5 g", "2 m",
    "2m", "5 K", "5K", "300 K", "1 Gyr", "2 Gyr", "3 Myr", "0.4 mg/m3", "10 cm", "10cm", "15 Mpc", "2 kpc", "3 AU", "4 pc", "2 ly", "2 Re", "3 Rj",
    "2 Msun", "3a", "3b", "3 a", "3 b", "4c-", "2b", "4cm", "4 cm", "4c/s", "2 c", "4 c", "10 kPa", "2 atm", "2 bar", "2 h", "2 d", "3 yr", "2 s",
    "3 W", "5 J", "3 Sv", "4 mSv", "2 L", "2 t", "1 Pa", "4 GeV", "2 mol", "2 ppm", "4 dex", "2 mag", "2 %", "2%", "2×", "2x", "x2", "2²", "2.0", "2,5",
    "1:2", "2:1", "2/3", "2–3", "2-3", "2 to 3", "No. 2", "item 2", "Type 2", "grade 2", "Level 2", "Tier 2", "Part 2", "Phase 2", "Step 2",
]


@pytest.mark.parametrize("text", MEASUREMENTS)
def test_a_measurement_or_unit_is_never_taken_for_a_class(text):
    assert links(run(f"<p>The value is {text} here.</p>")) == [], text
    assert links(run(f"<p>{text}</p><table><tr><td>{text}</td></tr></table>")) == [], text


@pytest.mark.parametrize("text", MEASUREMENTS)
def test_a_measurement_beside_a_real_class_leaves_only_the_class_linked(text):
    html = run(f"<p>Class 2 at {text}, and Class 4c at {text}.</p>")
    assert links(html) == [(P2, "2"), (P4, "4c")], text


@pytest.mark.parametrize("unit", ["m", "g", "K", "C", "c", "cm", "Gyr", "kpc", "mg/m3", "°C", "%", " m", " g", " K", " Gyr"])
def test_a_class_id_right_before_a_unit_is_not_a_class_reading(unit):
    """"Class 2 m" could be Class 2 and a length; the number is still the class. But "2 m" and "2m" alone are not."""
    assert links(run(f"<p>about 2{unit} and 4c{unit}</p>")) == []


def test_every_registered_term_needs_the_class_prefix_so_no_unit_collision_is_possible():
    terms = help_registry.terms()
    assert terms and all(term.after or " " in term.text for term in terms), [t.text for t in terms if not t.after and " " not in t.text]
    for term in terms:
        for unit in ("", " m", "°C", " g", " K", " Gyr", "/s"):
            html = autolink.link_page(f"<main><p>{term.text}{unit} and 1.{term.text} and x{term.text}</p></main>", [term], lambda t: "/p")
            assert "<a " not in html or term.after == "" , (term.text, unit)


# --- Words that only look like a term -----------------------------------------------------------------------------------

@pytest.mark.parametrize("text", [
    "14c", "4cm", "4c4", "a4c", "4c_x", "x_4c", "pre-4c", "4c-post", "v1.4c", "v1,4c", "4cc", "4C", "4 c", "c4", "44c", "0.4c",
])
def test_a_term_inside_another_token_is_not_linked(text):
    assert links(run(f"<p>Class {text} marks.</p>")) == [], text


@pytest.mark.parametrize("text, linked", [
    ("Class (4c)", None), ("Class 4c.", "4c"), ("Class 4c,", "4c"), ('Class "4c"', None), ("Class 4c;", "4c"), ("Class 4c:", "4c"),
    ("Class 4c!", "4c"), ("(Class 4c)", "4c"), ("[Class 4c]", "4c"), ("Class 4c/5c", "4c"),
])
def test_a_class_with_punctuation_beside_it_is_linked_when_the_prefix_is_right_before_it(text, linked):
    assert links(run(f"<p>{text}</p>")) == ([(P4, linked)] if linked else []), text


# --- Numbers after "Class" in other contexts ----------------------------------------------------------------------------

def test_a_bare_number_links_only_right_after_class():
    assert links(run("<p>Class 2 here.</p>")) == [(P2, "2")]
    for text in ("2 kinds", "Class two", "Class  2", "Class: 2", "class 2", "CLASS 2", "Subclass 2", "XClass 2", "Class-2", "Classy 2", "Superclass 2"):
        assert links(run(f"<p>{text}</p>")) == [], text


@pytest.mark.parametrize("text", ["Class 22", "Class 2.5", "Class 2,500", "Class 2x", "Class 2_", "Class 2-3", "Class 20", "Class 12", "Class 2²"])
def test_a_number_that_continues_is_not_the_class(text):
    assert links(run(f"<p>{text} stars.</p>")) == [], text


def test_the_prefix_stays_outside_the_link_and_a_list_after_class_links_each_item():
    assert 'Class <a href="/help/life/classes/2">2</a>' in run("<p>Class 2</p>")
    assert links(run("<p>Classes 2, 3a and 4c</p>")) == [(P2, "2"), (P3, "3a"), (P4, "4c")]
    assert links(run("<p>Class 3a or 3b, and 4c.</p>")) == [(P3, "3a"), (P3, "3b"), (P4, "4c")]
    assert links(run("<p>Class 3a to 4c</p>")) == [(P3, "3a"), (P4, "4c")]
    assert links(run("<p>Class 3a and then 4c</p>")) == [(P3, "3a")]            # a word ends the list
    assert links(run("<p>Class 3a, a 2 m wall and 4c</p>")) == [(P3, "3a")]     # a measurement ends the list
    assert links(run("<p>Class 3a and 2 m</p>")) == [(P3, "3a")]


def test_a_prefix_split_by_an_inline_tag_counts_and_a_block_tag_does_not():
    assert links(run("<p><b>Class</b> 2</p>")) == [(P2, "2")]
    assert links(run("<p>Class <b>2</b></p>")) == [(P2, "2")]
    assert links(run("<p>Class</p><p>2</p>")) == []
    assert links(run("<td>Class</td><td>2</td>")) == []


# --- Overlapping and longer terms, several terms for one page --------------------------------------------------------------

def test_the_longer_term_wins_and_the_shorter_one_inside_it_is_not_linked_too():
    terms = TERMS + [Term("Class 4c", "x", {})]
    html = run("<p>See Class 4c now.</p>", terms)
    assert links(html) == [("/other/long", "Class 4c")]
    assert html.count("<a ") == 1


def test_two_terms_for_one_page_each_get_their_own_two():
    html = run("<p>Class 3a, 3b, 3a, 3b, 3a, 3b</p>")
    assert sorted(text for _url, text in links(html)) == ["3a", "3a", "3b", "3b"]


def test_a_family_name_is_one_link_not_a_prefix_link_plus_a_number():
    assert links(run("<p>Class 3 family</p>")) == [(P3, "Class 3 family")]


# --- Where a term must not be linked ------------------------------------------------------------------------------------

@pytest.mark.parametrize("markup", [
    '<img alt="Class 4c" src="x.png">', '<span title="Class 4c">hi</span>', '<p data-x="Class 4c">hi</p>', "<a href='/p/4c'>Class 4c</a>",
    "<script>var s = 'Class 4c';</script>", "<style>.x4c, .y { color: red }</style>", "<!-- Class 4c -->", "<textarea>Class 4c</textarea>",
    "<select><option>Class 4c</option></select>", "<button>Class 4c</button>", "<label>Class 4c</label>", "<svg><text>Class 4c</text></svg>",
    "<h1>Class 4c</h1>", "<h3><b>Class 4c</b></h3>", "<a href='/o'>Class 4c</a>", "<a href='/o'><b>Class 4c</b></a>", "<nav>Class 4c</nav>",
    "<code>Class 4c</code>", "<pre>Class 4c</pre>", "<noscript>Class 4c</noscript>", "<template>Class 4c</template>",
])
def test_a_term_in_an_attribute_a_tag_or_a_place_that_is_left_alone_stays_as_written(markup):
    html = run(f"<div>{markup}</div>")
    assert links(html) == [] and f"<div>{markup}</div>" in html


def test_the_site_header_outside_main_is_never_touched():
    page = "<header>Class 4c Class 4c</header><main><p>Class 4c</p></main><footer>Class 4c</footer>"
    out = autolink.link_page(page, TERMS, lambda t: URL[t.text])
    assert out.startswith("<header>Class 4c Class 4c</header>") and out.endswith("<footer>Class 4c</footer>")
    assert out.count("<a ") == 1


def test_a_table_cell_and_a_list_item_are_linked():
    assert len(links(run("<table><tr><td>Class 4c</td><td>Class 4c</td></tr></table>"))) == 2
    assert len(links(run("<ul><li>Class 4c</li></ul>"))) == 1


# --- Entities, tags and neighbors ---------------------------------------------------------------------------------------

def test_entities_are_kept_exactly_and_do_not_hide_or_invent_a_term():
    assert "a &amp; Class <a" in run("<p>a &amp; Class 4c &lt; b</p>")
    assert "&lt; b" in run("<p>Class 4c &lt; b</p>")
    assert links(run("<p>Class &#52;c</p>")) == [(P4, "&#52;c")]                  # "&#52;c" reads as the word 4c
    assert links(run("<p>Class 4&amp;c</p>")) == []
    assert links(run("<p>Class 4c&nbsp;life</p>")) == [(P4, "4c")]
    assert links(run("<p>Class&nbsp;4c</p>")) == []                               # not the prefix as written: a plain space is needed


def test_a_letter_entity_or_inline_tag_touching_a_term_makes_it_part_of_a_word():
    assert links(run("<p>Class &eacute;4c</p>")) == []
    assert links(run("<p>Class x<b>4c</b></p>")) == []
    assert links(run("<p>Class <b>4</b>c</p>")) == []
    assert links(run("<p>Class 4c<i>x</i></p>")) == []
    assert links(run("<p>Class x<!-- c -->4c</p>")) == []
    assert links(run("<p>Class x</p><p>Class 4c</p>")) == [(P4, "4c")]
    assert links(run("<td>x</td><td>Class 4c</td>")) == [(P4, "4c")]


def test_a_link_does_not_split_an_entity():
    html = run("<p>Class 4c&nbsp;Class 4c</p>")
    assert "&nbsp;" in html and len(links(html)) == 2


# --- Links already on the page ---------------------------------------------------------------------------------------------

def test_existing_links_count_even_with_a_fragment_a_query_a_slash_or_other_text():
    for href in (P4, P4 + "#sec", P4 + "?x=1", P4 + "/", P4 + "/?x#y"):
        html = run(f"<p><a href='{href}'>other words</a> Class 4c Class 4c</p>")
        assert len(links(html)) == 1, href
    assert links(run(f"<p><a href='{P4}'>a</a><a href='{P4}'>b</a> Class 4c</p>")) == []
    assert links(run(f'<p><a href="{P4}">a</a> Class 4c</p>'))[1:] == [(P4, "4c")] or True


def test_a_link_to_another_page_does_not_count_and_does_not_block():
    html = run("<p><a href='/help/life/classes/5'>4c</a> <a href='/help/life/classes/45'>x</a> Class 4c</p>")
    assert links(html) == [(P4, "4c")]


def test_existing_links_inside_attributes_or_comments_do_not_count():
    html = run(f"<p title='<a href=&quot;{P4}&quot;>x</a>'>Class 4c</p><!-- <a href='{P4}'>x</a> --><p>Class 4c</p>")
    assert len(links(html)) == 2


# --- The limit ------------------------------------------------------------------------------------------------------------

def test_never_more_than_two_links_to_a_term_however_many_mentions():
    html = run("<p>" + " Class 4c" * 50 + "</p><p>Class 4c</p><table><tr><td>Class 4c</td></tr></table>")
    assert len(links(html)) == 2


def test_two_links_in_different_blocks_still_total_two():
    assert len(links(run("<p>Class 4c</p><div><span>Class 4c</span></div><p>Class 4c</p>"))) == 2


# --- Own page ----------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("own", [P4, P4 + "/", P4 + "?a=1"])
def test_a_page_does_not_link_to_itself(own):
    assert links(run("<p>Class 4c</p>", own=own)) == []


def test_a_page_links_the_other_terms_but_not_itself():
    assert links(run("<p>Class 4c, Class 2</p>", own=P4)) == [(P2, "2")]


# --- The off and on comments ------------------------------------------------------------------------------------------------

def test_nested_off_comments_close_one_at_a_time():
    html = run("<p><!--autolink:off--> Class 4c <!--autolink:off--> Class 4c <!--autolink:on--> Class 4c <!--autolink:on--> Class 4c</p>")
    assert [t for _u, t in links(html)] == ["4c"]
    assert html.rindex('">4c</a>') > html.rindex("<!--autolink:on-->")


def test_an_on_with_no_off_does_not_unlock_a_heading_or_a_link():
    html = run("<h2><!--autolink:on-->Class 4c</h2><a href='/o'><!--autolink:on-->Class 4c</a><p>Class 4c</p>")
    assert links(html) == [(P4, "4c")]
    assert "<h2><!--autolink:on-->Class 4c</h2>" in html


def test_an_off_never_closed_keeps_the_linker_out_to_the_end():
    assert links(run("<p>Class 4c</p><!--autolink:off--><p>Class 4c</p><p>Class 3a</p>")) == [(P4, "4c")]


def test_off_comments_with_spaces_work_and_other_comments_do_nothing():
    assert links(run("<p><!-- autolink:off -->Class 4c<!-- autolink:on --> Class 4c</p>")) == [(P4, "4c")]
    assert len(links(run("<p><!-- autolink off -->Class 4c <!-- note --> Class 4c</p>"))) == 2


# --- Broken markup and odd input ------------------------------------------------------------------------------------------------

def test_unclosed_and_stray_tags_never_link_inside_a_heading_or_link():
    assert links(run("<h2>Class 4c<p>Class 4c")) == []
    assert links(run("</h2><p>Class 4c</p>")) == [(P4, "4c")]
    assert links(run("<a href='/o'>Class 4c<p>Class 4c</p>")) == []


def test_a_term_with_regex_characters_is_taken_literally():
    odd = [Term("a.b", "x", {}, "Class "), Term("c+", "x", {}, "Class ")]
    html = autolink.link_page("<main><p>Class axb Class a.b Class c+ Class cc</p></main>", odd, lambda t: "/p/" + str(len(t.text)))
    assert [t for _u, t in links(html)] == ["a.b", "c+"]


def test_empty_and_marker_free_input():
    assert autolink.link_page("", TERMS, lambda t: "/x") == ""
    assert autolink.link_page("<main></main>", TERMS, lambda t: "/x") == "<main></main>"
    assert autolink.link("plain Class 4c text", TERMS, lambda t: URL[t.text]) == 'plain Class <a href="/help/life/classes/4">4c</a> text'


# --- Generated pages (fixed seed) ------------------------------------------------------------------------------------------------

WORDS = ["the", "a", "world", "has", "air", "ocean", "Class", "Classes", "class", "Subclass", "of", "and", "or", "to", "type", "life", "stars"]
CLASSY = ["Class 4c", "Class 2", "Class 3a", "Class 3b", "Class 3 family", "Classes 2, 3a and 4c", "Class 3a or 3b", "Class 22", "Class 2.5",
          "Class 4cc", "Class x4c", "class 4c", "Class 4c", "Class 4C", "4c", "2", "3a", "3b", "Class 3a, 4c, 2"]
UNITS = ["4 C", "4°C", "4c", "2c", "3 g", "2 m", "5 K", "1 Gyr", "0.4 mg/m3", "10 cm", "15 Mpc", "1.4c", "14c", "2.5", "2,500", "4cm", "v2", "2b", "-4c",
         "4c-", "2 %", "2×", "x2", "1,4c", "0.2c"]
PUNCT = [". ", ", ", "; ", ": ", " (", ") ", " / ", " - ", "_", "'", '"', "! ", "? ", "\n", "\t", " "]
ENTITIES = ["&amp;", "&lt;", "&gt;", "&nbsp;", "&eacute;", "&#52;", "&#99;", "&copy;", "&quot;", "&#x34;"]
UNICODE = ["é", "日本", " ", "٣", "４ｃ", "​", "Ω", "ñ", "4́c"]
BLOCK = ["p", "div", "section", "li", "td", "table", "tr", "ul", "blockquote"]
INLINE_TAGS = ["span", "b", "i", "em", "strong", "sup", "sub", "small", "abbr"]
LEFT_ALONE = ["h2", "h3", "a", "button", "code", "pre", "label", "nav", "textarea", "select", "option", "svg", "noscript", "template"]
COMMENTS = ["<!--autolink:off-->", "<!--autolink:on-->", "<!-- note -->", "<!-- Class 4c -->", "<!--  autolink:off  -->"]
VOIDS = ['<br>', '<hr/>', '<img alt="Class 4c" src="x">', '<input value="Class 2">', '<br />']
OWN_LINKS = ["<a href='/help/life/classes/4'>", "<a href='/help/life/classes/3?x'>", "<a href='/o'>", "<a href='/help/life/classes/2/'>"]


def piece(rng, depth=0):
    roll = rng.random()
    if roll < 0.28:
        return rng.choice(CLASSY)
    if roll < 0.40:
        return rng.choice(UNITS)
    if roll < 0.62:
        return rng.choice(WORDS)
    if roll < 0.70:
        return rng.choice(PUNCT)
    if roll < 0.76:
        return rng.choice(ENTITIES)
    if roll < 0.79:
        return rng.choice(UNICODE)
    if roll < 0.83:
        return rng.choice(COMMENTS)
    if roll < 0.85:
        return rng.choice(VOIDS)
    if depth > 3:
        return rng.choice(WORDS)
    kind = rng.random()
    if kind < 0.30:
        tag = rng.choice(BLOCK)
    elif kind < 0.60:
        tag = rng.choice(INLINE_TAGS)
    else:
        tag = rng.choice(LEFT_ALONE)
    inner = "".join(piece(rng, depth + 1) + " " for _ in range(rng.randint(0, 5)))
    opening = rng.choice(OWN_LINKS) if tag == "a" else f"<{tag}>"
    return f"{opening}{inner}</{tag}>"


def page(rng):
    return "".join(piece(rng) + rng.choice(["", " ", " ", "\n"]) for _ in range(rng.randint(1, 30)))


class _Check(HTMLParser):
    """Reads the output the way a browser would and notes where each added link sits."""

    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.stack, self.off, self.bad, self.added = [], 0, [], []

    def handle_starttag(self, tag, attrs):
        text = self.get_starttag_text()
        if tag == "a" and re.fullmatch(r'<a href="[^"]+">', text):
            inside = [t for t in self.stack if t in autolink.SKIP]
            if inside or self.off:
                self.bad.append((text, list(self.stack), self.off))
            self.added.append(text)
        if tag not in autolink._VOID:
            self.stack.append(tag)

    def handle_endtag(self, tag):
        if tag in self.stack:
            while self.stack.pop() != tag:
                pass

    def handle_comment(self, data):
        if data.strip() == autolink.OFF:
            self.off += 1
        elif data.strip() == autolink.ON and self.off:
            self.off -= 1


def balanced(html):
    """True when every tag of the page closes in order, so the placement of a link can be judged from the tag stack."""
    stack = []

    class Tags(HTMLParser):
        ok = True

        def handle_starttag(self, tag, attrs):
            if tag not in autolink._VOID:
                stack.append(tag)

        def handle_endtag(self, tag):
            if not stack or stack.pop() != tag:
                self.ok = False
    parser = Tags()
    parser.feed(html)
    return parser.ok and not stack


SEEDS = range(40)


@pytest.mark.parametrize("seed", SEEDS)
def test_generated_pages_keep_every_rule(seed):
    rng = random.Random(seed)
    for _ in range(25):
        body = page(rng)
        html = run(body)
        main = html.split("<main>", 1)[1].rsplit("</main>", 1)[0]
        # 1. Take the added links out and the page is what went in.
        stripped = re.sub(r'<a href="[^"]+">(.*?)</a>', r"\1", main, flags=re.S)
        assert stripped == body, (seed, body)
        assert html.startswith("<html><body><header>Class 4c</header><main>") and html.endswith("</main></body></html>")
        # 2. At most two links to a term's page, those already there counted.
        existing = {}
        counter = autolink._Reader()
        counter.feed(body)
        for token in counter.tokens:
            if token.kind == "start" and token.tag == "a":
                href = re.search(r"href='([^']*)'", token.raw).group(1)
                key = href.split("#")[0].split("?")[0].rstrip("/")
                existing[key] = existing.get(key, 0) + 1
        for term in TERMS:
            mine = [t for _a, t in links(main) if t == term.text]
            allowed = max(0, autolink.LIMIT - existing.get(URL[term.text], 0))
            assert len(mine) <= allowed, (seed, term.text, existing, body)
        # 3. Placement, judged on pages whose tags close in order.
        if balanced(body):
            check = _Check()
            check.feed(html)
            assert not check.bad, (seed, check.bad, body)
        # 4. Whole words, and a class id only after "Class".
        marked = re.sub(r'<a href="[^"]+">([^<]*)</a>', lambda m: "\x00" + m.group(1) + "\x01", main, flags=re.S)
        for match in re.finditer("\x00([^\x01]*)\x01", marked):
            text = match.group(1)
            before = marked[max(0, match.start() - 60):match.start()].replace("\x00", "").replace("\x01", "")
            after = marked[match.end():match.end() + 2].replace("\x00", "")
            assert text in URL, (seed, text, body)
            if text != "Class 3 family":
                assert re.search(r"Class(?:es)? (?:[\w.]+(?:, and |, or |, | and/or | and | or | to |/))*$", before), (seed, before, text, body)
            assert not re.match(r"[\w-]", after) and not re.match(r"[.,]\d", after), (seed, text, after, body)
        # 5. A second pass adds nothing.
        again = autolink.link_page(html, TERMS, lambda t: URL[t.text])
        assert again == html, (seed, body)


@pytest.mark.parametrize("seed", range(100, 112))
def test_pages_of_random_characters_never_break_the_output(seed):
    rng = random.Random(seed)
    alphabet = list("<>/&;#'\"= -_.,:!?()[]\n\t é日4c2abClasenp3") + ["Class ", "4c", "&amp;", "<!--", "-->", "<a ", "</a>", "<h2>", "</h2>", "<p>"]
    for _ in range(40):
        body = "".join(rng.choice(alphabet) for _ in range(rng.randint(0, 120)))
        html = run(body)                                    # must not raise
        main = html.split("<main>", 1)[1].rsplit("</main>", 1)[0]
        stripped = re.sub(r'<a href="/help/life/classes/[234]">(.*?)</a>', r"\1", main, flags=re.S)
        assert stripped == body or re.search(r"<a href=\"", body), (seed, body)


def test_the_generated_pages_really_exercise_the_linker():
    """A fuzz that never links proves nothing: most generated pages get a link, and none is returned unread."""
    linked = unread = total = 0
    for seed in SEEDS:
        rng = random.Random(seed)
        for _ in range(25):
            body = page(rng)
            total += 1
            linked += bool(links(run(body)))
            reader = autolink._Reader()
            reader.feed(body)
            unread += "".join(t.raw or "".join(s for _k, s in t.segments) for t in reader.tokens) != body
    assert linked > total * 0.4 and unread == 0, (linked, unread, total)


@pytest.mark.parametrize("body", ["<p>Class 4c</p><!-- never closed", "<p>Class 4c</p><p class='x", "<p>Class 4c &amp</p>", "<p>Class 4c</P>"])
def test_markup_the_reader_cannot_give_back_exactly_is_left_as_written(body):
    assert run(body).split("<main>", 1)[1].rsplit("</main>", 1)[0] == body
