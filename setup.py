# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='structure_factors',
    description='calculate the structure factors for grating interferometry',
    long_description=long_description,
    url="https://github.com/Enucatl/structure-factors",
    author='Spyros Gkoumas, Pablo Villanueva Perez and Matteo Abis',
    author_email='matteo.abis@psi.ch',
    license="MIT",
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
    ],
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=[
        'numpy',
        'scipy',
        'nist_lookup',
    ],
    entry_points="""
    [console_scripts]
    lynch_prediction = bin.lynch_prediction:main
    saxs_prediction = bin.saxs_prediction:main
    quick_plot = bin.quick_plot:main
    """
)
