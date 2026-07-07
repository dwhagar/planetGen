# stellarObjects/config.py

"""
Configuration settings for the stellarObjects package.

This module defines the `SystemConfig` class, which encapsulates all
configuration flags and variables that control the behavior of the star system
generation process. These settings can be modified via command-line arguments
in `planetGen.py` or directly within the code for specific scenarios.
"""

class SystemConfig:
    def __init__(self):
        self.system_flavor_count = 0
        self.recent_flavor_texts = [] # Stores recently used flavor texts to avoid repetition
        self.MARKDOWN = False
        """
        bool: If True, output will be formatted in Markdown. If False, output will be
        formatted in Wikitext. Defaults to False.
        """

        self.FORCE_HABITABLE_WORLD = False
        """
        bool: If True, the system generation will attempt to force the creation of at
        least one habitable world. Defaults to False.
        """

        self.FORCE_ASTEROID_BELT = False
        """
        bool: If True, the system generation will attempt to force the creation of at
        least one asteroid belt. Defaults to False.
        """

        self.FORCE_LARGE_STAR = False
        """
        bool: If True, the system generation will force the creation of a larger, more
        massive star. Defaults to False.
        """

        self.FORCE_MOONS = False
        """
        bool: If True, planets generated will have an increased likelihood of having
        moons. Defaults to False.
        """

        self.FORCE_MAX_PLANETS = False
        """
        bool: If True, the system generation will attempt to create the maximum number
        of planets the system can support. Defaults to False.
        """

        self.ABSURD = False
        """
        bool: If True, the system will generate an extremely large and dense system
        with the largest possible star, maximum planets, and moons. Defaults to False.
        """

        self.NO_PLANETS = False
        """
        bool: If True, the system generation will skip the planet generation process
        entirely, creating only a star. Defaults to False.
        """

        self.STAR_TYPE = None
        """
        str or None: Specifies a precise star type to generate (e.g., 'G2V'). If None,
        a star type will be randomly generated. Defaults to None.
        """

        self.NAME = None
        """
        str or None: Specifies a name for the star system, overriding the default
        random generation. If None, a name will be randomly generated. Defaults to None.
        """

        self.AGE = None
        """
        str or None: Specifies the age of the star system, either "young" or "old".
        If None, the age will be randomly determined. Defaults to None.
        """

        self.FORCE_INT = False
        """
        bool: If True, the system generation will attempt to force the creation of at
        least one planet with intelligent life. Implies FORCE_HABITABLE_WORLD.
        Defaults to False.
        """

        self.NO_INT = False
        """
        bool: If True, the system generation will ensure no planet with intelligent
        life is generated. Implies FORCE_HABITABLE_WORLD. Defaults to False.
        """

        self.NO_HABITABLE_WORLD = False
        """
        bool: If True, the system generation will ensure no habitable world is generated
        in the system. Defaults to False.
        """

        self.IS_BINARY_SYSTEM = False
        """
        bool: If True, the system generation will attempt to create a binary star system.
        Defaults to False.
        """