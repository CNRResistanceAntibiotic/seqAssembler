#!/usr/bin/python3

# setup.py file install libraries wanted for seqassembler
# launch it with pip install -e .

# Install setuptools if not already present.
from setuptools import setup, find_packages
import glob


def version():
    return "1.0.2"


setup(
    name='seqassembler',
    version=version(),
    description='seqassembler: pipeline CNR Resistance for assembling genomes',
    packages=find_packages(),
    author='Richard Bonnet',
    author_email='rbonnet@chu-clermontferrand.fr',
    python_requires='>=3.5',
    url='https://github.com/CNRResistanceAntibiotic/seqAssembler',
    scripts=glob.glob('scripts/*'),
    install_requires=[''],
    license='GPLv3',
    classifiers=[
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],

)
