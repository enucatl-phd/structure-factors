import numpy as np


def d_prime_yashiro(d_prime):
    """
    Function to calculate the geometrical factor of the DFEC
    Coefficient according to
    Appl Opt. 2011 Aug 1; 50(22): 4310–4319
    for d_prime > 1
    :d_prime: Diameter of the spheres / correlation_length fo the GI
    """
    return (
        d_prime -
        np.sqrt(d_prime ** 2 - 1) * (1 + 1 / (2 * d_prime ** 2)) +
        (1 / d_prime - 1 / (4 * d_prime ** 3)) *
        np.log(
            (d_prime + np.sqrt(d_prime ** 2 - 1)) /
            (d_prime - np.sqrt(d_prime ** 2 - 1))
        )
    )


def dark_field_extinction_coefficient(
        wavelength,
        grating_pitch,
        intergrating_distance,
        diameters,
        volume_fraction,
        delta_chi_squared):
    """
    Function to calculate the DFEC Coefficient
    For a spherical model from Appl Opt. 2011 Aug 1; 50(22):
    4310–4319

    :wavelength:
    :grating_pitch:
    :intergrating_distance:
    :diameters: diameter of the spheres, or vector of diameters
    :volume_fraction: fraction of the volume occupied by the microspheres
    :delta_chi_squared: difference in refractive index squared
    :returns: the dark field extinction coefficient μ_d

    """

    autocorrelation_length = wavelength * intergrating_distance / grating_pitch
    d_prime = diameters / autocorrelation_length
    factor = (3 * np.pi ** 2 * volume_fraction *
              delta_chi_squared *
              autocorrelation_length / wavelength ** 2)
    try:
        d_prime[d_prime > 1] = d_prime_yashiro(d_prime[d_prime > 1])
    except TypeError:
        # raised if d_prime is a float and not an array
        if d_prime > 1:
            d_prime = d_prime_yashiro(d_prime)
    return factor * d_prime
