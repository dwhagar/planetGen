# stellarObjects/names.py

"""
Name Generation Constants
=========================

This module contains all the data required for the procedural name generation
of celestial bodies, including base names, prefixes, suffixes, and word lists
for validation.
"""

import os
import nltk
from nltk.corpus import words

# --- Phonetic and Syllable Constants ---
VOWELS = "aeiou"
BAD_CONSONANTS = ["sz", "dt", "bp", "fb", "td", "pb"]

# --- Base Name Lists ---

# A list of names to be used as a base for generating star names.
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

# A list of names to be used as a base for generating planet names.
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

# A list of names to be used as a base for generating moon names.
MOON_NAMES = [
    "Dryad", "Naiad", "Nereid", "Oread", "Calliope", "Clio", "Erato", "Euterpe", "Melpomene", 
    "Polyhymnia", "Terpsichore", "Thalia", "Urania", "Charon", "Iris", "Janus", "Pan", "Brynhildr", 
    "Gunnr", "Hildr", "Sigrun", "Urd", "Verdandi", "Skuld", "Ratatoskr", "Vedrfolnir", "Aine", 
    "Banshee", "Puca", "Sidhe", "Amatsu", "Kunitsu", "Suijin", "Tenjin", "Kodama", "Tengu", "Kappa", 
    "Urvashi", "Menaka", "Rambha", "Tilottama", "Nandi", "Bhringi", "Chandesha", "Kokopelli", "Wendigo"
]

# --- Word Dictionaries and Validation ---

# Load the NLTK words corpus
nltk.download('words', quiet=True)
DICTIONARY_WORDS = set(words.words())

# Calculate the mean word size for splitting long words
WORD_SIZE_MEAN = round(sum(len(word) for word in DICTIONARY_WORDS) / len(DICTIONARY_WORDS))

# Load a list of NSFW words from a file.
with open(os.path.join(os.path.dirname(__file__), 'offensive_words.txt'), 'r') as f:
    NSFW_WORDS = {line.strip() for line in f}

# A list of common English words to avoid generating as names.
AVOIDED_NAMES = {
    "and", "are", "for", "you", "not", "the", "all", "new", "was", "can", "has", "but", "our", "one", "may", 
    "out", "use", "any", "see", "his", "who", "now", "get", "how", "its", "top", "had", "day", "two", "buy", 
    "her", "add", "jan", "she", "set", "map", "way", "off", "did", "car", "own", "end", "him", "per", "big", 
    "law", "art", "usa", "old", "non", "why", "low", "man", "job", "too", "men", "box", "air", "yes", "hot", 
    "say", "dec", "san", "tax", "got", "let", "act", "red", "key", "few", "age", "oct", "pay", "war", "nov", 
    "fax", "yet", "sun", "run", "net", "put", "try", "god", "log", "faq", "fun", "sep", "lot", "ask", "due", 
    "mar", "pro", "aug", "ago", "apr", "via", "bad", "far", "jun", "oil", "from", "that", "this", "with", 
    "your", "have", "more", "will", "home", "page", "free", "time", "they", "site", "what", "news", "only", 
    "when", "here", "also", "help", "view", "been", "were", "some", "like", "than", "find", "date", "back", 
    "list", "name", "just", "over", "year", "into", "next", "used", "work", "last", "most", "data", "make", 
    "them", "post", "city", "such", "best", "then", "good", "info", "high", "each", "very", "book", "read", 
    "need", "many", "user", "said", "does", "mail", "full", "life", "know", "days", "part", "real", "item", 
    "must", "made", "line", "send", "type", "take", "area", "want", "long", "code", "show", "even", "much", 
    "sign", "file", "link", "open", "case", "same", "both", "game", "care", "down", "size", "shop", "text", 
    "rate", "form", "love", "john", "main", "time", "year", "people", "way", "day", "man", "thing", "woman", 
    "life", "child", "world", "school", "state", "family", "student", "group", "country", "problem", "hand", 
    "part", "place", "case", "week", "company", "system", "program", "question", "work", "government", 
    "number", "night", "point", "home", "water", "room", "mother", "area", "money", "story", "fact", 
    "month", "lot", "right", "study", "book", "eye", "job", "word", "business", "issue", "side", "kind", 
    "head", "house", "service", "friend", "father", "power", "hour", "game", "line", "end", "member", 
    "law", "car", "city", "community", "name", "president", "team", "minute", "idea", "kid", "body", 
    "information", "back", "parent", "face", "others", "level", "office", "door", "health", "person", 
    "art", "war", "history", "party", "result", "change", "morning", "reason", "research", "girl", "guy", 
    "moment", "air", "teacher", "force", "education"
}

# --- Name Affixes ---

# Prefixes for star names
STAR_PREFIXES = ["Al", "El", "Il", "Ul", "O", "E", "A", "I"]
# Suffixes for star names
STAR_SUFFIXES = ["ia", "a", "os", "us", "is", "es", "e", "o"]

# Prefixes for planet names
PLANET_PREFIXES = ["Ze", "Xe", "Ve", "Ge", "Pe", "Te", "Ke", "Re"]
# Suffixes for planet names
PLANET_SUFFIXES = ["a", "i", "o", "u", "ia", "io", "iu", "ea"]

# Prefixes for moon names
MOON_PREFIXES = ["Li", "Mi", "Ni", "Pi", "Si", "Ti", "Ki", "Ri"]
# Suffixes for moon names
MOON_SUFFIXES = ["a", "e", "i", "o", "u"]
