# stellarObjects/asteroidData.py

"""
Asteroid Belt Generation
========================

This module contains the `AsteroidBelt` class, which provides a representation
for asteroid belts within a star system. It also includes utility functions
related to planet mass ranges, which are used by both planets and asteroid belts.
"""

import random
from . import config

# Define common asteroid components
ASTEROID_COMPONENTS = [
    "carbon", "oxygen", "silicon", "magnesium", "aluminum", "calcium",
    "sulfur", "phosphorus", "hydrogen", "nitrogen", "iron", "nickel",
    "iridium", "palladium", "platinum", "gold", "osmium", "ruthenium",
    "rhodium", "olivine", "pyroxene", "plagioplase feldspars", "kamacite",
    "taenite", "troilite", "schreibersite", "cohenite", "serpentine",
    "magnetite", "hematite", "chromite", "carbon monoxide", "carbon dioxide",
    "silicon carbide", "iron sulfide"
]


class AsteroidBelt:
    """
    A basic class to store information for an asteroid belt.

    This class serves as a simple container for the properties of an asteroid
    belt, primarily its orbital distance and boundaries.

    Attributes:
        distance (float): The average distance of the asteroid belt from the star in AU.
        lower_limit (float): The inner boundary of the asteroid belt in AU.
        upper_limit (float): The outer boundary of the asteroid belt in AU.
        type (str): A character representing the object type, 'a' for asteroid belt.
        density (str): The density of the asteroid belt ('dense', 'sparse', 'typical').
        composition (list): A list of common compounds found in the belt.
    """

    def __init__(self, distance, lower_limit, upper_limit):
        """
        Initializes an AsteroidBelt object.

        Args:
            distance (float): The average distance from the star in AU.
            lower_limit (float): The inner boundary of the belt in AU.
            upper_limit (float): The outer boundary of the belt in AU.
        """
        self.distance = distance
        self.lower_limit = lower_limit
        self.upper_limit = upper_limit
        self.type = 'a'
        self.density = random.choice(["dense", "sparse", "typical"])
        self.composition = self._generate_composition()

    def _generate_composition(self):
        """
        Generates a list of common compounds for the asteroid belt.
        Includes at least 3 random components, plus frozen water and hydrocarbons.
        """
        selected_components = random.sample(ASTEROID_COMPONENTS, k=min(3, len(ASTEROID_COMPONENTS)))

        # Ensure frozen water and hydrocarbons are always included
        if "frozen water" not in selected_components:
            selected_components.append("frozen water")
        if "various long and short chain hydrocarbons" not in selected_components:
            selected_components.append("various long and short chain hydrocarbons")

        return selected_components

    def to_paragraph_list(self):
        """
        Returns a list of strings, where each string is a paragraph describing
        the asteroid belt.
        """
        # Conversion factor from AU to Light-Years. 1 LY = 63241.1 AU
        AU_TO_LY = 1 / 63241.1
        LY_THRESHOLD = 1.0  # The distance in LY at which to switch from AU to LY display

        # Convert the belt's boundaries to light-years to check against the threshold
        upper_limit_ly = self.upper_limit * AU_TO_LY

        if upper_limit_ly < LY_THRESHOLD:
            # For smaller systems, display the boundaries in AU for better precision
            distance_text = f"between {self.lower_limit:.3f} AU and {self.upper_limit:.3f} AU"
        else:
            # For very large systems, display in light-years for readability
            lower_limit_ly = self.lower_limit * AU_TO_LY
            distance_text = f"between {lower_limit_ly:.4f} light-years and {upper_limit_ly:.4f} light-years"

        output_paragraphs = []
        header_level = '##' if config.MARKDOWN else '=='
        header = f"{header_level} Asteroid Belt {header_level if not config.MARKDOWN else ''}"
        output_paragraphs.append(header)

        composition_text = ", ".join(self.composition)
        description_sentence = f"This is a {self.density} asteroid belt, primarily composed of {composition_text}."

        output_paragraphs.append(f"An asteroid belt orbits roughly {distance_text}. {description_sentence}")

        return output_paragraphs

    def __str__(self):
        """
        Returns a string representation of the asteroid belt, with paragraphs
        separated by double newlines.
        """
        return "\n\n".join(self.to_paragraph_list())