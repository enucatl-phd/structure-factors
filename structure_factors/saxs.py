import numpy as np


def sphere_form_factor(q, d):
    """
    Function to calculate the sphere form factor for different radii hard
    spheres
    D.J. Kinning et al., Macromolecules 17 (1984) 1712

    :q: 4*pi/lambda*sin(theta/2)
    :d: diameter of the microspheres
    """
    a = q * d / 2
    form_factor = 3 / a ** 3 * (np.sin(a) - a * np.cos(a))
    form_factor[a == 0] = 1 / 3
    return form_factor


def hard_sphere_structure_factor(q, d, f):
    """
    Function to calculate structure factors for different radii hard spheres
    D.J. Kinning et al., Macromolecules 17 (1984) 1712

    :q: 4*pi/lambda*sin(theta/2)
    :d: diameter of the microspheres
    :f: fraction volume
    """

    a = q * d

    alpha = ((1 + 2 * f) ** 2) / ((1 - f) ** 4)
    beta = -6 * f * (1 + (f / 2)) ** 2 / (1 - f) ** 4
    gamma = 0.5 * f * (1 + (2 * f) ** 2) / (1 - f) ** 4

    G1 = alpha / (a ** 2) * (np.sin(a) - (a * np.cos(a)))
    G2 = beta / (a ** 3) * ((2 * a * np.sin(a)) +
                            ((2 - (a ** 2)) * np.cos(a)) - 2)
    G3 = gamma / (a ** 5) * (((-a ** 4) * np.cos(a)) +
                             (4 * (((3 * a ** 2 - 6) * np.cos(a)) +
                                   (np.sin(a) * (a ** 3 - 6 * a)) + 6)))
    G = (G1 + G2 + G3)

    S = 1 / (1 + (24 * f * (G / a)))
    S[a == 0] = 1 / (1 + 24 * f * (alpha / 3 + beta / 4 + gamma / 6))
    return S


def dark_field_extinction_coefficient(
        wavelength,
        grating_pitch,
        intergrating_distance,
        diameter,
        volume_fraction,
        delta_chi_squared,
        real_space_sampling):
    """TODO: Docstring for autocorrelation.
    # n = even

    :arg1: TODO
    :returns: TODO

    """

    autocorrelation_length = wavelength * intergrating_distance / grating_pitch
    sampling_step = real_space_sampling[1] - real_space_sampling[0]
    n = len(real_space_sampling)
    fourier_sampling_step = 1 / (n * sampling_step)
    fourier_space_sampling = np.linspace(
        -n / (16 * autocorrelation_length),
        n / (16 * autocorrelation_length),
        n,
        endpoint=False,
    )
    fx, fy = np.meshgrid(fourier_space_sampling, fourier_space_sampling)
    f = np.sqrt(fx ** 2 + fy ** 2)
    q = 2 * np.pi * f

    sphere_factor = sphere_form_factor(q, diameter)
    structure_factor = hard_sphere_structure_factor(
        q, diameter, volume_fraction
    )

    intensity = (
        volume_fraction *
        delta_chi_squared *
        np.pi * diameter ** 3 / 6 *
        sphere_factor * structure_factor
    )

    k = 2 * np.pi / wavelength
    autocorrelation = np.real(
        k ** 2 *
        np.fft.ifft2(np.fft.ifftshift(intensity)) *
        fourier_sampling_step ** 2
    )

    autocorrelation_sampling = np.fft.ifftshift(real_space_sampling)
    autocorrelation_values = autocorrelation[0, :]

    autocorrelation1 = np.interp(
        autocorrelation_length,
        autocorrelation_sampling,
        autocorrelation_values
    )

    return autocorrelation_values[0] - autocorrelation1
