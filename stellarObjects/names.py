# stellarObjects/names.py

"""
Name Generation Constants
=========================

This module contains all the data required for the procedural name generation
of celestial bodies, including base names, prefixes, suffixes, and word lists
for validation. The data is organized to provide a rich and varied source for
creating unique and plausible-sounding names for stars, planets, moons, and
sectors.
"""

import os

import nltk
from nltk.corpus import words

# --- Phonetic and Syllable Constants ---

VOWELS = "aeiou"
"""A string containing the vowels of the English alphabet, used for syllable splitting."""

BAD_CONSONANTS = ["sz", "dt", "bp", "fb", "td", "pb", "sz", "zs", "srs", "srl", "szl", "szr", "zsl", "zsr", "sch", "rsc", "lsc", "rsh", "lsh", "tsh", "dsh", "ksh", "gsh", "csh", "xsh", "qsh", "psh", "bsh", "fsh", "vsh", "msh", "nsh", "ssh", "zsh", "hsh", "wsh", "jsh", "ysch", "rsch", "lsch", "tsch", "dsch", "ksch", "gsch", "csch", "xsch", "qsch", "psch", "bsch", "fsch", "vsch", "msch", "nsch", "ssch", "zsch", "hsch", "wsch", "jsch"]
"""
A list of consonant clusters that are considered phonetically unpleasing or
difficult to pronounce in generated names. This list helps ensure that the
names are pronounceable and sound more natural.
"""

# --- Base Name Lists ---

STAR_NAMES = [
    "Mimosa", "Gacrux", "El Nath", "Miaplacidus", "Alnilam", "Alnair", "Alioth", "Dubhe", "Mirfak", 
    "Wezen", "Sargas", "Kaus", "Australis", "Avior", "Menkalinan", "Alkaid", "Ras", "Algethi", "Alhena", 
    "Peacock", "Alsuhail", "Navi", "Regor", "Alphard", "Algieba", "Gienah", "Markab", "Enif", 
    "Scheat", "Alpheratz", "Mirach", "Diphda", "Hamal", "Sheratan", "Mesarthim", "Caph", "Schedar", 
    "Ruchbah", "Marfak", "Albireo", "Izar", "Kochab", "Pherkad", "Thuban", "Eltanin", "Rastaban", 
    "Grumium", "Kuma", "Zosma", "Chara", "Cor", "Caroli", "Vindemiatrix", "Porrima", "Auva", "Heze", 
    "Zavijava", "Syrma", "Unukalhai", "Zubenelgenubi", "Zubeneschamali", "Shaula", "Lesath", "Nunki", 
    "Kaus", "Borealis", "Kaus", "Media", "Alnasl", "Ascella", "Facies", "Rukbat", "Arkab", "Prior", 
    "Arkab", "Posterior", "Algedi", "Dabih", "Nashira", "Sadalmelik", "Sadalsuud", "Skat", "Ancha", 
    "Situla", "Formalhaut", "Hadar", "Rigil", "Kentaurus", "Toliman", "Agena", "Proxima", 
    "Centauri", "Alpha", "Centauri", "Beta", "Centauri", "Barnard's", "Wolf", "Lalande", "Epsilon", 
    "Eridani", "Lacaille", "Luyten", "Tau", "Ceti", "Gliese", "Kapteyn's", "Kepler", "Helios", "Aten", 
    "Ra", "Horus", "Osiris", "Isis", "Thoth", "Anubis", "Seth", "Hathor", "Bastet", "Sekhmet", "Ptah", 
    "Amun", "Khonsu", "Sobek", "Taweret", "Maat", "Nut", "Geb", "Shu", "Tefnut", "Zeus", "Hera", "Poseidon", 
    "Hades", "Demeter", "Hestia", "Apollo", "Artemis", "Ares", "Aphrodite", "Hephaestus", "Hermes", 
    "Dionysus", "Athena", "Persephone", "Eros", "Nike", "Morpheus", "Hypnos", "Thanatos", "Nemesis", 
    "Tyche", "Iris", "Hecate", "Juno", "Phoebus", 
    "Diana", "Vulcan", "Bacchus", "Minerva", "Proserpina", "Cupid", "Victoria", 
    "Somnus", "Mors", "Invidia", "Fortuna", "Arcus", "Hecuba", "Odin", "Frigg", "Thor", "Loki", "Baldur", 
    "Heimdall", "Tyr", "Freya", "Freyr", "Njord", "Skadi", "Idunn", "Bragi", "Sif", "Hel", "Fenrir", 
    "Jormungandr", "Surtr", "Ymir", "Dagda", "Lugh", "Brigid", "Manannan", "Morrigan", "Aengus", "Dian", 
    "Cecht", "Nuada", "Goibniu", "Midir", "Danu", "Epona", "Cernunnos", "Taranis", "Belenus", "Esus", 
    "Teutates", "Raiden", "Fujin", "Susanoo", "Amaterasu", "Tsukuyomi", "Izanagi", "Izanami", "Hachiman", 
    "Inari", "Bishamon", "Benten", "Daikoku", "Ebisu", "Fudo", "Kannon", "Jizo", "Marici", "Yasha", 
    "Ashura", "Brahma", "Vishnu", "Shiva", "Saraswati", "Lakshmi", "Parvati", "Ganesha", "Kartikeya", 
    "Indra", "Agni", "Surya", "Chandra", "Varuna", "Vayu", "Yama", "Kubera", "Kali", "Durga", "Hanuman", 
    "Rama", "Krishna", "Sita", "Radha", "Zodiac", "Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo", 
    "Libra", "Scorpio", "Sagittarius", "Capricorn", "Aquarius", "Pisces",
    "Anu", "Enlil", "Enki", "Inanna", "Marduk", "Tiamat", "Apsu", "Kingu", "Ereshkigal", "Nergal",
    "Huitzilopochtli", "Quetzalcoatl", "Tezcatlipoca", "Tlaloc", "Tonatiuh", "Metztli", "Coyolxauhqui", "Xipe", "Totec", "Mictlantecuhtli",
    "Itzamna", "Kukulkan", "Hunab", "Ku", "Ixchel", "Ah", "Puch", "Kinich", "Ahau", "Chaac", "Yum", "Kaax",
    "Viracocha", "Inti", "Mama", "Quilla", "Pachamama", "Supay", "Illapa", "Kon", "Catequil", "Apu",
    "Olorun", "Olodumare", "Obatala", "Shango", "Yemoja", "Ogun", "Oshun", "Orunmila", "Eshu", "Aganju",
    "Chukwu", "Ala", "Amadioha", "Ekwensu", "Ikenga", "Anyanwu", "Agbara", "Idemili", "Ogbunabali", "Ojukwu",
    "Bondye", "Legba", "Kalfu", "Marinet", "Damballa", "Ayida", "Wedo", "Erzulie", "Freda", "Ogou", "Feray",
    "Bunjil", "Gnowee", "Baiame", "Daramulum", "Yhi", "Lumaluma", "Marmoo", "Ngintaka", "Tjinimin", "Wuriupranili",
    "Rangi", "Papa", "Tangaroa", "Tane", "Mahuta", "Maui", "Hina", "Pele", "Lono", "Ku",
    "Ukko", "Akka", "Tapio", "Mielikki", "Ahti", "Vellamo", "Ilmarinen", "Louhi", "Pekko", "Lemminkainen"
]
"""
A list of names to be used as a base for generating star names. These names are
drawn from various mythologies, real star names, and other astronomical sources
to provide a diverse and epic feel.
"""

PLANET_NAMES = [
    "Nix", "Hydra", "Styx", "Kerberos", "Adrastea", 
    "Amalthea", "Thebe", "Metis", "Aitne", "Ananke", "Aoede", "Arche", "Autonoe", "Callirrhoe", "Carme", 
    "Carpo", "Chaldene", "Cyllene", "Dia", "Eirene", "Elara", "Erinome", "Ersa", "Euanthe", "Eukelade", 
    "Eupheme", "Euporie", "Eurydome", "Harpalyke", "Hegemone", "Helike", "Hermippe", "Herse", "Himalia", 
    "Iocaste", "Isonoe", "Kale", "Kallichore", "Kalyke", "Kore", "Leda", "Lysithea", "Megaclite", 
    "Melpomene", "Mneme", "Orthosie", "Pandia", "Pasiphae", "Pasithee", "Philophrosyne", "Praxidike", 
    "Sinope", "Sponde", "Taygete", "Thelxinoe", "Themisto", "Thyone", "Valetudo", "Aegaeon", "Aegir", 
    "Albiorix", "Alcyone", "Atlas", "Bebhionn", "Beli", "Bergelmir", "Bestla", "Calypso", "Daphnis", 
    "Epimetheus", "Erriapus", "Farbauti", "Fenrir", "Fornjot", "Greip", "Hati", "Helene", "Hyperion", 
    "Hyrrokkin", "Ijiraq", "Janus", "Jarnsaxa", "Kari", "Kiviuq", "Loge", "Methone", "Mundilfari", 
    "Narvi", "Paaliaq", "Pallene", "Pan", "Pandora", "Phoebe", "Polydeuces", "Prometheus", "Siarnaq", 
    "Skathi", "Skoll", "Surtur", "Suttungr", "Tarqeq", "Tarvos", "Telesto", "Thrymr", "Ymir", 
    "Belinda", "Bianca", "Caliban", "Cordelia", "Cressida", "Cupid", "Desdemona", "Ferdinand", "Francisco", 
    "Juliet", "Mab", "Margaret", "Ophelia", "Perdita", "Portia", "Prospero", "Puck", "Rosalind", "Setebos", 
    "Stephano", "Sycorax", "Trinculo", "Despina", "Galatea", "Halimede", "Hippocamp", "Laomedeia", "Larissa", 
    "Naiad", "Neso", "Proteus", "Psamathe", "Sao", "Thalassa", "Namaka", "Hi'iaka",
    "Cahokia", "Chaco", "Mesa", "Verde", "Taos", "Acoma", "Hovenweep", "Kinishba", "Paquime", "Snaketown",
    "Cherokee", "Apache", "Navajo", "Sioux", "Iroquois", "Hopi", "Zuni", "Seminole", "Choctaw", "Cree",
    "Gitche", "Manitou", "Wakan", "Tanka", "Coyote", "Raven", "Iktomi", "Glooscap", "Nanabozho", "Hiawatha",
    "Teotihuacan", "Tikal", "Copan", "Palenque", "Tenochtitlan", "Cusco", "Machu", "Picchu", "Nazca", "Caral",
    "Aztec", "Maya", "Inca", "Olmec", "Zapotec", "Mixtec", "Toltec", "Taino", "Arawak", "Carib",
    "Quetzalcoatl", "Itzamna", "Viracocha", "Huitzilopochtli", "Hunab", "Ku", "Inti", "Kukulkan", "Tlaloc", "Tezcatlipoca"
]
"""
A list of names to be used as a base for generating planet names. This list
includes names of moons from our solar system, as well as names from various
mythologies and ancient cultures, providing a wide range of naming conventions.
"""

MOON_NAMES = [
    "Dryad", "Naiad", "Nereid", "Oread", "Calliope", "Clio", "Erato", "Euterpe", "Melpomene", 
    "Polyhymnia", "Terpsichore", "Thalia", "Urania", "Charon", "Iris", "Janus", "Pan", "Brynhildr", 
    "Gunnr", "Hildr", "Sigrun", "Urd", "Verdandi", "Skuld", "Ratatoskr", "Vedrfolnir", "Aine", 
    "Banshee", "Puca", "Sidhe", "Amatsu", "Kunitsu", "Suijin", "Tenjin", "Kodama", "Tengu", "Kappa", 
    "Urvashi", "Menaka", "Rambha", "Tilottama", "Nandi", "Bhringi", "Chandesha", "Kokopelli", "Wendigo"
]
"""
A list of names to be used as a base for generating moon names. These names are
often mythological or whimsical, suitable for smaller celestial bodies orbiting
a planet.
"""

# --- Word Dictionaries and Validation ---

# Load the NLTK words corpus
nltk.download('words', quiet=True)
DICTIONARY_WORDS = set(words.words())
"""
A set of common English words from the NLTK corpus. This is used to validate
generated names and ensure they are not actual words, which helps in creating
unique and fictional-sounding names.
"""

WORD_SIZE_MEAN = round(sum(len(word) for word in DICTIONARY_WORDS) / len(DICTIONARY_WORDS))
"""
The calculated mean (average) length of words in the `DICTIONARY_WORDS` set.
This value is used to determine if a generated name is long enough to be split
into two parts for better readability.
"""

with open(os.path.join(os.path.dirname(__file__), 'offensive_words.txt'), 'r') as f:
    NSFW_WORDS = {line.strip() for line in f}
"""
A set of "Not Safe For Work" (NSFW) or offensive words, loaded from an external
file. This list is used as a filter to prevent the generation of inappropriate
or offensive names.
"""

# --- Name Affixes ---

STAR_PREFIXES = [
    "Al", "El", "Il", "Ul", "O", "E", "A", "I",  # Existing
]
"""
A list of prefixes for star names. These prefixes are designed to sound grand,
important, or celestial. They are drawn from various languages to add a sense
of diversity and scale to the star names.
"""

STAR_SUFFIXES = [
    "ia", "a", "os", "us", "is", "es", "e", "o",  # Existing
    "um", "or", "ex", "ion", "ius", "ae", "oris", "ax",  # Latin
    "dras", "cyon", "nar", "tor"  # Other
]
"""
A list of suffixes for star names. Similar to the prefixes, these suffixes
contribute to a grand and epic-sounding name, suitable for stars.
"""

PLANET_PREFIXES = [
    "Ze", "Xe", "Ve", "Ge", "Pe", "Te", "Ke", "Re",  # Existing
    "Verd", "Azul", "Rojo", "Blanc", "Noir",  # Colors
    "Erd", "Stein", "Wasser", "Wald",  # German
    "Piedra", "Agua", "Bosque"  # Spanish/Italian
]
"""
A list of prefixes for planet names. These prefixes are generally themed around
terrestrial concepts like colors and natural elements (earth, stone, water),
providing a more grounded feel compared to star names.
"""

PLANET_SUFFIXES = [
    "a", "i", "o", "u", "ia", "io", "iu", "ea",  # Existing
    "ara", "eth", "os", "ica", "or", "es", "ana", "is",  # Common
    "dine", "tine", "don", "gan"  # Other
]
"""
A list of suffixes for planet names. These suffixes often give a sense of place
or commonality, fitting for planetary bodies.
"""

MOON_PREFIXES = [
    "Li", "Mi", "Ni", "Pi", "Si", "Ti", "Ki", "Ri",  # Existing
    "Whis", "Ech", "Sha", "Glim", "Gleam", "Fae", "Pix", "Spri", "Wyn",  # Whimsical
    "Luni", "Umb", "Somni", "Ani", "Magi", "Rune",  # Mystical
    "Traum", "Geist", "Seel",  # German
    "Rêv", "Ombr",  # French
    "Somb", "Sueñ", "Sogn"  # Spanish/Italian
]
"""
A list of prefixes for moon names. The theme for these prefixes is whimsical and
mystical, reflecting the smaller, often more enchanting nature of moons.
"""

MOON_SUFFIXES = [
    "a", "e", "i", "o", "u",  # Existing
    "elle", "ette", "ina", "ia", "is", "ie",  # Diminutive/Feminine
    "ix", "yx", "ax", "ex", "ra", "la", "sa"  # Mystical
]
"""
A list of suffixes for moon names. These suffixes complement the mystical and
whimsical prefixes, often with diminutive or feminine endings to create
delicate-sounding names.
"""

SECTOR_NAMES = [
    # Galactic arms, superclusters, and other real large-scale structures
    "Perseus", "Sagittarius", "Orion", "Norma", "Scutum", "Centaurus", "Cygnus",
    "Virgo", "Laniakea", "Coma", "Saraswati", "Shapley", "Hercules", "Fornax",
    "Eridanus", "Ophiuchus", "Hydra", "Pisces", "Cetus", "Gould",
    # Real nebulae and star clusters
    "Helix", "Tarantula", "Trifid", "Rosette", "Cocoon", "Lagoon", "Omega",
    "Veil", "Boomerang", "Flame", "Horsehead", "Carina", "Pleiades", "Hyades",
    "Trapezium", "Beehive", "Jewel", "Butterfly",
    # Warhammer 40,000 sectors and expanses
    "Calixis", "Koronus", "Jericho", "Scarus", "Ixaniad", "Obscurus", "Pacificus",
    "Tempestus",
    # EVE Online regions
    "Deklein", "Delve", "Derelik", "Detorid", "Devoid", "Esoteria", "Etherium",
    "Everyshore", "Exordium", "Fade", "Feythabolis", "Fountain", "Geminate",
    "Genesis", "Heimatar", "Immensea", "Impass", "Insmother", "Kador", "Kalevala",
    "Khanid", "Lonetrek", "Malpais", "Metropolis", "Oasa", "Omist", "Paragon",
    "Perrigen", "Placid", "Pochven", "Providence", "Querious", "Solitude",
    "Stain", "Syndicate", "Tenal", "Tenerifis", "Tribute", "Venal", "Yasna",
    # Star Wars regions and sectors
    "Corellia", "Kuat", "Hutt", "Tion", "Chommell", "Bpfassh",
    # Star Trek regions
    "Briar", "Badlands", "Mutara", "Typhon", "Denorios",
    # Asimov's Foundation series sectors
    "Arcturus", "Normannic", "Trantor", "Terminus", "Anacreon", "Korell",
    "Periphery",
    # Mass Effect regions
    "Attican", "Argos", "Traverse",
    # Elite Dangerous regions
    "Formidine", "Odin",
    # Dune
    "Arrakis", "Corrino", "Atreides", "Harkonnen", "Kaitain", "Landsraad",
    # Firefly and Babylon 5
    "Verse", "Georgia", "Narn", "Minbar"
]
"""
A list of names to be used as a base for generating sector names. These are
drawn from real astronomical large-scale structures (galactic arms,
superclusters, nebulae, and star clusters) alongside space-region and sector
names from science fiction (Warhammer 40,000, EVE Online, Star Wars, Star
Trek, Asimov's Foundation, Mass Effect, Elite Dangerous, Dune, Firefly, and
Babylon 5), giving generated sectors the same epic, well-traveled feel as the
star and planet name lists.
"""

SECTOR_PREFIXES = [
    "Kal", "Kor", "Sca", "Zar", "Vel", "Dra", "Xan", "Nor",  # Existing-style short affixes
    "Tar", "Vor", "Kry", "Bel", "Fen", "Rhu", "Il", "Os"  # Frontier/expanse-flavored
]
"""
A list of prefixes for sector names. These are short, hard-edged affixes
meant to evoke the frontier-sector naming style of works like Warhammer
40,000 and EVE Online.
"""

SECTOR_SUFFIXES = [
    "is", "us", "on", "ar", "ine", "eth",  # Existing-style short affixes
    "yx", "oth", "ade", "ium", "ain", "ost"  # Expanse/reach-flavored
]
"""
A list of suffixes for sector names, giving the generated names the same
grand, place-name feel as the base `SECTOR_NAMES` they're combined with.
"""

UNIVERSAL_PHONEMES = [
    # Romance (Spanish/Italian/French/Portuguese)
    "rr", "ll", "nh", "lh", "gn", "ci", "gi", "eu", "oi", "ch",
    # Germanic (German/Dutch/Scandinavian)
    "pf", "ij", "aa", "oo", "ae", "oe", "ck", "tz",
    # Slavic (Russian/Polish/Czech, transliterated)
    "zh", "shch", "ts", "cz", "rz", "ya", "yu", "yo", "ye",
    # Arabic/Hebrew/Persian (transliterated, no diacritics)
    "kh", "gh", "dh", "th", "qa", "ee", "ain",
    # Turkish
    "yz", "ogh",
    # South Asian (Hindi/Sanskrit/Bengali, transliterated)
    "bh", "dh", "gh", "kh", "ph", "chh", "jh", "jn", "shr",
    # East Asian - Mandarin (pinyin)
    "zh", "sh", "xi", "qi", "ang", "ong", "iong",
    # Japanese (romaji)
    "tsu", "kya", "gyo", "nya", "ryo", "hyu", "fu",
    # Korean (romanization)
    "eo", "eu", "ae", "oe", "ng", "gyu", "hye",
    # Vietnamese
    "ng", "nh", "tr", "gia", "qu",
    # Austronesian (Tagalog/Indonesian/Malay)
    "ngg", "ny", "mba",
    # Polynesian (Maori/Hawaiian/Samoan)
    "wh", "nga", "ae", "oa", "ai",
    # African (Swahili/Bantu)
    "mb", "nd", "nj", "gw", "vwa",
    # Mesoamerican (Nahuatl/Mayan)
    "tl", "tza", "huo", "xo", "cua", "ix", "b'a",
    # Andean (Quechua/Aymara, glottalized)
    "ch'", "q'o", "k'i",
    # Celtic (Irish/Welsh/Scottish Gaelic)
    "bh", "dh", "mh", "wy", "ao",
    # Finno-Ugric (Finnish/Hungarian/Estonian)
    "gy", "ty", "cs", "nye",
    # Caucasian (Georgian/Armenian)
    "kh", "dz", "gho",
]
"""
A pool of short, ASCII-7-bit-printable phoneme chunks -- rough romanized
approximations of consonant/vowel combinations drawn from as many of the
world's language families as practical (Romance, Germanic, Slavic,
Arabic/Hebrew/Persian, Turkish, South Asian, Mandarin, Japanese, Korean,
Vietnamese, Austronesian, Polynesian, Bantu, Mesoamerican, Andean, Celtic,
Finno-Ugric, and Caucasian), each written using only plain ASCII letters
and apostrophes (no diacritics, since diacritics fall outside 7-bit ASCII).

`generate_phoneme_salad_name` (see `stellarObjects/utils.py`) occasionally
splices one of these into the shuffled-syllable pool for every kind of
name -- star, planet, moon, and sector alike, since they all funnel
through that one shared function -- widening the cultural "naming budget"
generated names draw from beyond whatever happens to appear in that
type's own base `_NAMES` list.
"""
