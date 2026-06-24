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
from . import constants

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
        composition (list): A list of tuples, each containing (component, concentration).
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
        Generates a list of common compounds for the asteroid belt,
        including random components with unique, ordered concentrations.
        Returns a list of (component, concentration) tuples.
        """
        all_concentrations = ["high", "moderate", "small", "trace"]
        
        # Determine how many components to select (up to 4)
        num_components_to_select = min(len(all_concentrations), len(constants.ASTEROID_COMPONENTS))
        
        # Select unique components
        selected_components = random.sample(constants.ASTEROID_COMPONENTS, k=num_components_to_select)
        
        # Shuffle selected components to randomize which concentration they get
        random.shuffle(selected_components)
        
        # Take only the necessary number of concentrations
        concentrations_for_use = all_concentrations[:num_components_to_select]
        
        composition_data = []
        for i in range(num_components_to_select):
            component = selected_components[i]
            concentration = concentrations_for_use[i]
            composition_data.append((component, concentration))

        return composition_data

    def to_paragraph_list(self):
        """
        Returns a list of strings, where each string is a paragraph describing
        the asteroid belt.
        """
        # Convert the belt's boundaries to light-years to check against the threshold
        upper_limit_ly = self.upper_limit * constants.AU_TO_LY

        if upper_limit_ly < constants.LY_THRESHOLD:
            # For smaller systems, display the boundaries in AU for better precision
            distance_text = f"between {self.lower_limit:.3f} AU and {self.upper_limit:.3f} AU"
        else:
            # For very large systems, display in light-years for readability
            lower_limit_ly = self.lower_limit * constants.AU_TO_LY
            distance_text = f"between {lower_limit_ly:.4f} light-years and {upper_limit_ly:.4f} light-years"

        output_paragraphs = []
        header_level = '##' if config.MARKDOWN else '=='
        header = f"{header_level} Asteroid Belt {header_level if not config.MARKDOWN else ''}"
        output_paragraphs.append(header)

        composition_phrases = []
        for component, concentration in self.composition:
            if concentration == "trace":
                composition_phrases.append(f"trace amounts of {component}")
            else:
                composition_phrases.append(f"{concentration} concentrations of {component}")

        composition_text = ""
        if composition_phrases:
            if len(composition_phrases) == 1:
                composition_text = composition_phrases[0]
            elif len(composition_phrases) == 2:
                composition_text = f"{composition_phrases[0]} and {composition_phrases[1]}"
            else:
                composition_text = ", ".join(composition_phrases[:-1]) + f", and {composition_phrases[-1]}"
        else:
            composition_text = "various unknown materials" # Fallback if no components are selected

        description_sentence = f"This is a {self.density} asteroid belt, composed of {composition_text}."

        output_paragraphs.append(f"An asteroid belt orbits roughly {distance_text}. {description_sentence}")

        return output_paragraphs

    def __str__(self):
        """
        Returns a string representation of the asteroid belt, with paragraphs
        separated by double newlines.
        """
        return "\n\n".join(self.to_paragraph_list())