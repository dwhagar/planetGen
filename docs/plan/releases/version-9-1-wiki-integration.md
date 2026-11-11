# Version 9.1: The game's manual, linked

Filed on 2026-11-08 as a version map (Boss, 7:03 PM EST: a wiki-style site, spakr, will hold the game's manual, characters, races and the other information that PlanetGen does not store; PlanetGen is the game's starmap; a later, unstarted project will link a body to its terrain and map).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 9.1

Status: waiting on the wiki site (Boss, 7:07 PM EST on 2026-11-08: hold 9.1 until spakr has its public API and its field search); nothing is started or owned in PlanetGen, and 8.1 comes first.

## Goal

A visitor on an object page can reach the wiki pages written about that star, planet, moon, belt, comet or phenomenon, and the wiki can link back to the map. PlanetGen stores no wiki content and the wiki stores no map data: each side keeps its own accounts (Boss, 7:08 PM EST: separate accounts; an optional account sync is the later item USR.32) and the two sites connect with API keys made by the owner of each side, each holding the other's API URL (Boss, 7:06 PM EST). Every link names a body by the galaxy key and the Object ID, so a wiped map with a new seed reads as "from another map" and a link to a rebuilt system reads as "no longer charted". The map and every page keep working when the wiki is off, slow or down. Nothing from the removed wiki code (OPS.69) comes back, the spakr repository is not edited from here (the request to its owner is DOC.119), and a later terrain and map project answers the same lookup (API.47). Design: the draft `design-drafts/wiki-integration-spakr-2026-11-08.md` in the project files, until Boss approves it into docs/design.

## Phase 1: Identity and access

The galaxy key, the short address and the read-only integration keys that a linked site needs.

### Galaxy key and addresses

| ID | Item | Needs | Note |
|---|---|---|---|
| API.43 | A galaxy key and the site identity call |  |  |
| API.44 | The short address /o/<Object ID> that finds the kind and redirects, and one stable answer from /api/objects/<ref> that tells 'no longer charted' from 'wrong map'; check the anonymous read rule | API.43 |  |

### Integration keys

| ID | Item | Needs | Note |
|---|---|---|---|
| SEC.47 | Integration keys | USR.28 |  |
| ADM.152 | An Owner-only Integration keys screen: make, name, show once, rotate and revoke a key for a linked site, with its last use and audit | SEC.47 |  |

## Phase 2: The settings and the lookup

The wiki settings section, the lookup client and the linked-pages call.

### Settings and lookup

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.153 | The wiki settings section in config.json |  |  |
| API.45 | The lookup client | ADM.153 |  |
| API.46 | GET /api/objects/<ref>/linked-pages: the wiki pages about an object, asked of the wiki lazily and returned short | API.45 |  |

## Phase 3: The pages

What a visitor sees on the map.

### Wiki panel

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.296 | The Wiki panel on object and phenomenon pages, the mark on system tree rows, and the Help menu and footer link | API.46 |  |

## Phase 4: Tests, developer page and the request to spakr

The checks that the link cannot harm the map, and the documents for the wiki site's owner.

### Tests

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.710 | Tests that the wiki link can never harm the map: wiki down or slow, wrong galaxy key, GM-only pages never listed, and no wiki text in MySQL | UX.296 |  |

### Documents

| ID | Item | Needs | Note |
|---|---|---|---|
| DOC.118 | A developer page for spakr's author | API.43, API.44 |  |
| DOC.119 | A request to the spakr owner: a typed planetgen_ref field on place pages, a links lookup for public pages with a batch form, and an API version header | DOC.118 |  |

### Linked sites list

| ID | Item | Needs | Note |
|---|---|---|---|
| API.47 | The linked-sites list for the third project | API.45 |  |
