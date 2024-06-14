#!/usr/bin/env python

from setuptools import setup, find_packages


def readfile(filename):
    with open(filename, 'r+') as f:
        return f.read()


setup(name='tilepy',
      version='2.0.0',
      description='Computation of the tiling scheduling of large localization uncertainty region event with multi-wavelength pointing telescopes',
      install_requires=[
          'astropy', #>=4.1',
          'scipy', #',
          'healpy', #==1.12.9',
          'ipython', #',
          'matplotlib<3.9.0', #==2.2.2',
          'MOCpy', #==0.10.0',
          'numpy', #==1.19.5',
          'pandas',
          'pytz', #==2021.1',
          'ephem', #==3.7.6.0',
          'gdpyc',
          'tables',
          'fastparquet',
          'skyfield',
          'ligo.skymap', #~=1.0.7',
      ],
      packages=find_packages(),
      # tests_require=['pytest'],
      author='Monica Seglar-Arroyo, Halim Ashkar, Fabian Schussler, Mathieu de Bony ',
      author_email='mosear2@gmail.com ',
      url='https://github.com/astro-transients/tilepy',
      long_description=readfile('README.md'),
      include_package_data=True,
      license='GNU3',
      classifiers=[],
      scripts=[],
      )
