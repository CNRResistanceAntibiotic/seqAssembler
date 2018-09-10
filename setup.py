#!/usr/bin/python3

# setup.py file install libraries wanted for seqassembler
# launch it with pip install -e .

# Install setuptools if not already present.
from setuptools import setup, find_packages
import glob

setup(
    name='seqassembler',
    version='1.0.1',
    description='seqassembler: pipeline CNR Resistance for assembl genomes',
    packages=find_packages(),
    author='Richard Bonnet',
    author_email='rbonnet@chu-clermontferrand.fr',
    url='https://github.com/CNRResistanceAntibiotic/seqAssembler',
    scripts=glob.glob('scripts/*'),
    install_requires=['pandas', 'matplotlib', 'pysam==0.11.2.2', 'pysamstats'],
    license='GPLv3',
    classifiers=[
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],

)
