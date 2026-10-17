# planetgen/physics/neutron_star.py

"""
The radius of a neutron star from its mass (GEN.463).

It replaces the fixed 10 to 13 km range, drawn uniformly with a floor the radius measurements do not support. The
radius of a neutron star is nearly independent of its mass between 1.4 and 2.1 solar masses, and what is known is a
band: Miller et al. 2021 (ApJL 918, L28, arXiv:2105.06979; PSR J0740+6620 with NICER and the gravitational-wave and
nuclear constraints) give R(1.4) = 12.45 +- 0.65 km and R(2.08) = 12.35 +- 0.75 km, each the full range of the 68%
intervals over three equation-of-state frameworks. The centre of the radius is taken on the straight line between those
two masses (carried on past them, as a nearly flat band), its spread likewise, and a star's radius is the centre plus a
normal draw of that spread, clipped to `tuning.NEUTRON_STAR_RADIUS_CLIP_SIGMA` standard deviations so that every
radius lies in the 95% band. The measurements it is checked against: PSR J0030+0451 at 1.34 solar masses, 12.71
+1.14 -1.19 km (Riley et al. 2019, arXiv:1912.05702), and PSR J0740+6620 at 2.08 solar masses, 13.7 +2.6 -1.5 km (Miller
et al. 2021).
"""

from planetgen import tuning


def radius_band_km(mass_solar):
    """The centre and spread, km, of the radius of a neutron star of this mass: the line through the two anchors of
    `tuning.NEUTRON_STAR_RADIUS_ANCHORS`."""
    (mass_low, centre_low, sigma_low), (mass_high, centre_high, sigma_high) = tuning.NEUTRON_STAR_RADIUS_ANCHORS
    fraction = (mass_solar - mass_low) / (mass_high - mass_low)
    return centre_low + (centre_high - centre_low) * fraction, sigma_low + (sigma_high - sigma_low) * fraction


def radius_km(mass_solar, deviation=0.0):
    """
    The radius, km, of a neutron star of this mass: the centre of the band plus `deviation` standard deviations of
    it, `deviation` (a standard normal draw; 0 for the centre) clipped to `NEUTRON_STAR_RADIUS_CLIP_SIGMA`.
    """
    limit = tuning.NEUTRON_STAR_RADIUS_CLIP_SIGMA
    centre, sigma = radius_band_km(mass_solar)
    return centre + sigma * min(max(deviation, -limit), limit)


def radius_limits_km(mass_solar):
    """The lowest and highest radius, km, `radius_km` gives for this mass."""
    return radius_km(mass_solar, -tuning.NEUTRON_STAR_RADIUS_CLIP_SIGMA), radius_km(mass_solar,
                                                                                    tuning.NEUTRON_STAR_RADIUS_CLIP_SIGMA)
