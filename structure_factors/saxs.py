import numpy as np
import scipy
import logging

log = logging.getLogger(__name__)


def sphere_form_factor(q, d):
    """
    Function to calculate the sphere form factor for different radii hard
    spheres
    D.J. Kinning et al., Macromolecules 17 (1984) 1712

    :q: 4 * pi / lambda * sin(theta / 2)
    :d: diameter of the microspheres
    """
    a = q * d / 2
    a = np.ma.array(a)
    form_factor = 3 / a ** 3 * (np.sin(a) - a * np.cos(a))
    return form_factor.filled(1 / 3)


def g(a, f):
    alpha = (1 + 2 * f) ** 2 / (1 - f) ** 4
    beta = -6 * f * (1 + (f / 2)) ** 2 / (1 - f) ** 4
    gamma = 0.5 * f * (1 + (2 * f)) ** 2 / (1 - f) ** 4
    a = np.ma.array(a)
    result = (
        alpha / (a ** 2) * (np.sin(a) - (a * np.cos(a))) +
        beta / (a ** 3) * (
            (2 * a * np.sin(a)) +
            ((2 - (a ** 2)) * np.cos(a)) - 2) +
        gamma / (a ** 5) * (
            -a ** 4 * np.cos(a) +
            (4 * (((3 * a ** 2 - 6) * np.cos(a)) +
                  ((a ** 3 - 6 * a) * np.sin(a)) + 6))
        )
    )
    result /= a
    return result.filled(alpha / 3 + beta / 4 + gamma / 6)


def hard_sphere_structure_factor(q, d, f):
    """
    Function to calculate structure factors for different radii hard spheres
    D.J. Kinning et al., Macromolecules 17 (1984) 1712

    :q: 4*pi/lambda*sin(theta/2)
    :d: diameter of the microspheres
    :f: fraction volume
    """
    sf = 1 / (1 + (24 * f * g(q * d, f)))
    log.debug("hard sphere structure factor \n %s", sf)
    return sf


def dark_field_extinction_coefficient(
        wavelength,
        grating_pitch,
        intergrating_distance,
        diameter,
        volume_fraction,
        delta_chi_squared,
        real_space_sampling,
        structure_factor_function=lambda *args: 1):
    """
    Function to calculate the DFEC Coefficient

    https://www.nature.com/articles/srep35259

    :wavelength:
    :grating_pitch:
    :intergrating_distance:
    :diameter: diameter of the spheres
    :volume_fraction: fraction of the volume occupied by the microspheres
    :delta_chi_squared: difference in refractive index squared
    :real_space_sampling: sampling in real space, recommended value
    np.linspace(
    -4 * autocorrelation_length,
    4 * autocorrelation_length,
    2048,
    endpoint=False,
    )
    :structure_factor_function: a function of q, diameter, volume fraction
    that calculates the relevant structure factor (defaults to 1)
    :returns: the dark field extinction coefficient μ_d

    """

    autocorrelation_length = wavelength * intergrating_distance / grating_pitch
    sampling_step = real_space_sampling[1] - real_space_sampling[0]
    log.debug("sampling step %g", sampling_step)
    n = len(real_space_sampling)
    log.debug("sampling cells %g", n)
    log.debug("real space sampling %s", real_space_sampling)
    fourier_sampling_step = 1 / (n * sampling_step)
    log.debug("fourier sampling step %g", fourier_sampling_step)
    fourier_space_sampling = np.linspace(
        -n / (16 * autocorrelation_length),
        n / (16 * autocorrelation_length),
        n,
        endpoint=False,
    )
    log.debug("fourier space sampling %s", fourier_space_sampling)
    fx, fy = np.meshgrid(fourier_space_sampling, fourier_space_sampling)
    f = np.sqrt(fx ** 2 + fy ** 2)
    q = 2 * np.pi * f

    sphere_factor = sphere_form_factor(q, diameter) ** 2
    log.debug("sphere form factor \n %s", sphere_factor)
    structure_factor = structure_factor_function(q, diameter, volume_fraction)

    intensity = (
        volume_fraction *
        delta_chi_squared *
        np.pi * diameter ** 3 / 6 *
        sphere_factor * structure_factor
    )
    log.debug("intensity \n %s", intensity)

    k = 2 * np.pi / wavelength
    log.debug("k \n %g", k)
    autocorrelation = (
        k ** 2 *
        np.fft.fftshift(np.real(np.fft.ifft2(np.fft.ifftshift(intensity)))) /
        sampling_step ** 2
    )
    log.debug("autocorrelation \n %s", autocorrelation)

    autocorrelation_sampling = np.fft.ifftshift(real_space_sampling)
    autocorrelation_values = np.fft.ifftshift(autocorrelation)[0, :]
    autocorrelation1 = scipy.interpolate.interp1d(
        autocorrelation_sampling,
        autocorrelation_values
    )(autocorrelation_length)
    log.debug("autocorrelation length \n %s", autocorrelation_length)
    log.debug("autocorrelation sampling \n %s", autocorrelation_sampling)
    log.debug("autocorrelation values \n %s", autocorrelation_values)
    log.debug("autocorrelation1 %g", autocorrelation1)

    dfec = autocorrelation_values[0] - autocorrelation1

    log.debug("saxs dfec %g", dfec)
    return dfec
