import numpy as np


def sphere_form_factor(q, r):
    """
    Function to calculate the sphere form factor for different radii hard spheres
    D.J. Kinning et al., Macromolecules 17 (1984) 1712

    :q: 4*pi/lambda*sin(theta/2)
    :r: radius of microspheres
    """
    a = q * r
    form_factor = 3 / a ** 3 * (np.sin(a) - a * np.cos(a))
    form_factor[a == 0] = 1 / 3
    return form_factor

def hard_sphere_structure_factor(q, r, f):
    """
    Function to calculate structure factors for different radii hard spheres
    D.J. Kinning et al., Macromolecules 17 (1984) 1712

    :q: 4*pi/lambda*sin(theta/2)
    :r: radius of microspheres
    :f: fraction volume
    """
    
    a = 2 * q * r

    alpha = ((1 + 2 * f) ** 2)/((1 - f) ** 4)
    beta = -6 * f * (1 + (f / 2)) ** 2 / (1 - f) ** 4
    gamma = 0.5 * f * (1 + (2 * f) ** 2) / (1 - f) ** 4

    G1 = alpha / (a ** 2) * (np.sin(a) - (a * np.cos(a)))
    G2 = beta / (a ** 3) * ((2 * a * np.sin(a)) + ((2 - (a ** 2)) * np.cos(a)) - 2)
    G3 = gamma / (a ** 5) * (((-a ** 4) * np.cos(a))
        + (4 * (((3 * a ** 2 - 6) * np.cos(a))
        + (np.sin(a) * (a ** 3 - 6 * a)) + 6)))
    G = (G1 + G2 + G3)

    S = 1 / (1 + (24 * f * (G/a)))
    S[a==0] = 1 / (1 + 24 * f * (alpha / 3 + beta / 4 + gamma / 6))
    return S
