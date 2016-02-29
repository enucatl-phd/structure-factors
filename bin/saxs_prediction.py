import click
import csv
import numpy as np
from scipy import constants

import structure_factors.saxs as saxs
from nist_lookup import xraydb_plugin as xdb


@click.command()
@click.option("--energy", type=float, default=25,
              help="design energy of the interferometer [keV]")
@click.option("--grating_pitch", type=float, default=2e-6,
              help="pitch of G2 [m]")
@click.option("--intergrating_distance", type=float, default=12e-2,
              help="pitch of G2 [m]")
@click.option("--min_diameter", type=float, default=0.1e-6,
              help="minimum diameter of the spheres [m]")
@click.option("--max_diameter", type=float, default=10e-6,
              help="maximum diameter of the spheres [m]")
@click.option("--diameter_step", type=float, default=0.1e-6,
              help="maximum diameter of the spheres [m]")
@click.option("--volume_fraction", type=float, default=0.4,
              help="fraction of the total volume occupied by the spheres")
@click.option("--sphere_material", default="SiO2",
              help="chemical composition of the spheres")
@click.option("--sphere_density", type=float, default=2.0,
              help="density of the material of the spheres [g/cm³]")
@click.option("--background_material", default="C3H8O3",
              help="chemical composition of the background")
@click.option("--background_density", type=float, default=1.26,
              help="density of the material of the background [g/cm³]")
@click.option("--output", type=click.File("w"), default="-",
              help="output file for the csv data")
def main(
        energy,
        grating_pitch,
        intergrating_distance,
        min_diameter,
        max_diameter,
        diameter_step,
        volume_fraction,
        sphere_material,
        sphere_density,
        background_material,
        background_density,
        output
        ):
    wavelength = (constants.physical_constants["Planck constant in eV s"][0] *
                  constants.c / (energy * 1e3))
    diameters = np.arange(min_diameter, max_diameter, diameter_step)
    delta_sphere, beta_sphere, _ = xdb.xray_delta_beta(
        sphere_material,
        sphere_density,
        energy * 1e3)
    delta_background, beta_background, _ = xdb.xray_delta_beta(
        background_material,
        background_density,
        energy * 1e3)
    delta_chi_squared = (
        (delta_sphere - delta_background) ** 2 +
        (beta_sphere - beta_background) ** 2
    )

    autocorrelation_length = wavelength * intergrating_distance / grating_pitch
    n = 2048
    real_space_sampling = np.linspace(
        -4 * autocorrelation_length,
        4 * autocorrelation_length,
        n,
        endpoint=False,
    )
    output_csv = csv.writer(output)
    output_csv.writerow(["diameter", "dfec"])
    for diameter in diameters:
        dfec = saxs.dark_field_extinction_coefficient(
            wavelength,
            grating_pitch,
            intergrating_distance,
            diameter,
            volume_fraction,
            delta_chi_squared,
            real_space_sampling
        )
        output_csv.writerow([diameter, dfec])
