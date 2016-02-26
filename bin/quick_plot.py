import click
import numpy as np
import matplotlib.pyplot as plt


@click.command()
@click.argument("inputfile", type=click.File("r"))
def main(inputfile):
    a = np.loadtxt(inputfile, delimiter=",", skiprows=1)
    plt.plot(a)
    plt.show(block=True)
