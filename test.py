from stellarObjects import Star, Planet

def main():
    """Creates a Sun-like star and an Earth-like planet for testing."""

    # Create a Sun-like star (G2V spectral class, 5778 K temperature)
    sun = Star(spectral_class="G", temperature=5778)

    # Calculate Earth's orbital distance in kilometers for the atmospheric pressure calculation
    earth_distance_km = 149.6e6  # 1 AU in kilometers

    # Create an Earth-like planet (radius 1 Earth radius, distance 1 AU)
    earth = Planet(sun.habitable_zone, 1, sun.luminosity, sun.radius, sun.temperature)

    earth.distance = 1 # in AU
    earth.planet_class = "M"
    earth.radius = 6371  # km
    earth.mass = 5.972e24  # kg
    earth.density = 5.515 # g/cm^3
    earth.atm_density = 1.225  # kg/m³
    earth.atm_molar_density = 0.02897 # kg/mol
    earth.atmosphere = "Nitrogen, Oxygen, Argon"
    earth.composition = "Silicon, Iron, Magnesium"

    earth.calculate_surface_gravity()
    earth.calculate_atmospheric_pressure()
    earth.calculate_atmospheric_temperature()

    # Calculate and print Earth's properties
    print(earth)
    print(f"Earth's Surface Gravity: {earth.gravity:.2f} g")
    print(f"Earth's Atmospheric Pressure: {earth.atmospheric_pressure:.2f} Pa")
    print(f"Earth's Surface Temperature: {earth.surface_temperature:.2f} K")

if __name__ == "__main__":
    main()