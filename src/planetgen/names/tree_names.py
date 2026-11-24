# planetgen/names/tree_names.py

"""
The names of the members, orbit nodes and planet hosts of a system's tree (GEN.312; extends GEN.62).

A tree is nodes (`orbit_nodes`: `node_no`, `parent_node_no`, `slot`, `kind`) and members hung in them (a star or a dark member: `(parent_node_no, slot)`).
The names are a pure function of that structure and the system's name (and, for a tree with a wide split, the branch words), so two builds of one seed give the
same names, and a stored system's names can be read back to the words they were made from.

- **A wide split** (a node of kind `wide`) gives each of its two children a branch word of its own, drawn from the seeded stream (`word_slots` says how many and in
  what order: the slot-0 child, the heavier, from the star word lists and the slot-1 child from the small/child ones, as GEN.62 does for two). The first word of the
  system's name is shared.
- **Inside a branch, a close split** gives each child a path code: level one `A`, `B`; level two `a`, `b`; level three `1`, `2`; deeper levels `.1`, `.2` appended
  (a `group` node's children go on past `B`: `C`, `D`). A star alone in its branch has no code, and a branch restarts the codes.
- **A star** is named `<base> <branch word> <code>` (the parts it has): `Voranthis Kelmoor Aa`, `Voranthis Pikkita Ba`, `Voranthis Ba` for a tree with no wide
  split. `base` is the system's name for a tree with no wide split and its first word for one with. A tree of one or two stars is named as GEN.62 names it: a close
  pair is `Voranthis A` and `Voranthis B`, a wide pair `Voranthis Kelmoor` and `Voranthis Pikkita`.
- **Planets** take a stem that is a word, never a bare letter: the branch that holds the primary uses the first word, any other branch its own word, and a coded
  star adds its code (`Voranthis I`, `Pikkita I`, `Pikkita Ba I`). The planets of a pair about its node (P-type) take the node's stem: the branch word at the
  branch root, else `<branch word> <code>`, and the system's name for the root of a tree with no wide split.
"""

from dataclasses import dataclass

LEVEL_ONE = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
LEVEL_TWO = "abcdefghijklmnopqrstuvwxyz"


def _get(item, key):
    return item[key] if isinstance(item, dict) else getattr(item, key)


def code_letter(level, index):
    """The code a child of index `index` (0 first) takes at path level `level` (1 first): `A`, `a`, `1`, then `.1` at the deeper levels. A letter level runs out at
    26 children; a digit level does not."""
    if level == 1:
        return LEVEL_ONE[index]
    if level == 2:
        return LEVEL_TWO[index]
    return str(index + 1) if level == 3 else f".{index + 1}"


@dataclass(frozen=True)
class TreeNames:
    """
    The names of one tree. `members[i]` is the name and `member_stems[i]` the planet stem of member `i` (in the order of the `members` argument),
    `nodes[node_no]` the name of the pair at a node and `node_stems[node_no]` the stem of the planets about it.
    """
    members: tuple
    member_stems: tuple
    nodes: dict
    node_stems: dict


def _children(nodes, members):
    """`{node_no: [(slot, "node", node_no) | (slot, "member", index), ...]}` in slot order."""
    found = {_get(node, "node_no"): [] for node in nodes}
    for node in nodes:
        parent = _get(node, "parent_node_no")
        if parent is not None:
            found[parent].append((_get(node, "slot"), "node", _get(node, "node_no")))
    for index, (parent, slot) in enumerate(members):
        found[parent].append((slot, "member", index))
    return {number: sorted(items) for number, items in found.items()}


def _root(nodes):
    roots = [_get(node, "node_no") for node in nodes if _get(node, "parent_node_no") is None]
    if len(roots) != 1:
        raise ValueError(f"a tree has one root node, this has {len(roots)}")
    return roots[0]


def word_slots(nodes, members):
    """
    The branch words the tree needs, as `[(node_no, slot), ...]` in the order they are drawn (nodes from the root, depth first, the slot-0 child first): one for each
    child of each wide node. A tree with no wide node needs none, a wide pair two (the primary's, then the secondary's).
    """
    if not nodes:
        return []
    kinds = {_get(node, "node_no"): _get(node, "kind") for node in nodes}
    children = _children(nodes, members)
    found = []

    def walk(number):
        for slot, kind, ref in children[number]:
            if kinds[number] == "wide":
                found.append((number, slot))
            if kind == "node":
                walk(ref)

    walk(_root(nodes))
    return found


def name_tree(system_name, nodes, members, words=(), primary=0):
    """
    The `TreeNames` of a tree. `nodes` are anything with `node_no`, `parent_node_no`, `slot` and `kind` (an `OrbitNode`, or a row as a dict); `members` is
    `[(parent_node_no, slot), ...]` for every star (`primary` is the index of the one whose branch uses the system's first word) and then every dark member; `words`
    are the branch words in the order `word_slots` gives (not needed for a tree with no wide node).
    """
    first_word = system_name.split()[0]
    if not nodes:
        return TreeNames((system_name,) * len(members), (system_name,) * len(members), {}, {})
    wide = word_slots(nodes, members)
    if len(words) != len(wide):
        raise ValueError(f"the tree needs {len(wide)} branch words, {len(words)} given")
    word_of = dict(zip(wide, words))
    kinds = {_get(node, "node_no"): _get(node, "kind") for node in nodes}
    children = _children(nodes, members)
    base = first_word if wide else system_name
    placed_members, placed_nodes = {}, {}

    def compose(*parts):
        return " ".join(part for part in parts if part)

    def walk(number, scope, word, code, level):
        for slot, kind, ref in children[number]:
            if kinds[number] == "wide":
                child = ((number, slot), word_of[(number, slot)], "")
                child_level = 1
            else:
                child = (scope, word, code + code_letter(level, slot))
                child_level = level + 1
            if kind == "member":
                placed_members[ref] = child
            else:
                placed_nodes[ref] = child
                walk(ref, *child, child_level)

    root = _root(nodes)
    walk(root, None, None, "", 1)
    primary_scope = placed_members[primary][0] if primary in placed_members else None

    def stem(scope, word, code):
        if word is None:
            return compose(base, code)
        return compose(first_word if scope == primary_scope else word, code)

    member_names = tuple(compose(base, *placed_members[index][1:]) for index in range(len(members)))
    member_stems = tuple(stem(*placed_members[index]) for index in range(len(members)))
    node_names = {root: system_name, **{number: compose(base, word, code) for number, (_scope, word, code) in placed_nodes.items()}}
    node_stems = {root: base, **{number: stem(*entry) for number, entry in placed_nodes.items()}}
    return TreeNames(member_names, member_stems, node_names, node_stems)


def words_from_names(system_name, nodes, members, names, primary=0, missing=None):
    """
    The branch words a stored tree's star names were made from, in `word_slots` order, or `None` where the names are not of that shape (a star renamed by hand,
    or a word that more than one name disagrees on). `names[i]` is the name of member `i`, `None` for a member with no name (a dark one). A branch no name is under
    (one that holds only dark members) has no word to read: `None`, or `missing` when it is given (a reader that names only nodes of lit stars).
    """
    wide = word_slots(nodes, members)
    if not wide:
        return ()
    placeholder = "\0"
    shape = name_tree(system_name, nodes, members, [placeholder] * len(wide), primary).members
    scopes = _scopes_of(nodes, members)
    found = {}
    for index, (stored, expected) in enumerate(zip(names, shape)):
        if stored is None or placeholder not in expected:
            continue
        head, _, tail = expected.partition(placeholder)
        word = stored[len(head):len(stored) - len(tail)] if stored.startswith(head) and stored.endswith(tail) else ""
        if not word or " " in word or len(head) + len(tail) >= len(stored):
            return None
        found[index] = word
    words = []
    for slot in wide:
        candidates = {found[index] for index, path in scopes.items() if slot in path and index in found}
        if not candidates and missing is not None:
            words.append(missing + str(len(words)))
        elif len(candidates) != 1:
            return None
        else:
            words.append(candidates.pop())
    return tuple(words) if len(set(words)) == len(words) else None


def _scopes_of(nodes, members):
    """`{member index: [(node_no, slot), ...]}`: every child of a wide node the member hangs under."""
    kinds = {_get(node, "node_no"): _get(node, "kind") for node in nodes}
    parent = {_get(node, "node_no"): (_get(node, "parent_node_no"), _get(node, "slot")) for node in nodes}
    found = {}
    for index, (at, slot) in enumerate(members):
        path = []
        while at is not None:
            if kinds[at] == "wide":
                path.append((at, slot))
            at, slot = parent[at]
        found[index] = path
    return found
