#!/usr/bin/env python

from setuptools import setup, find_packages


def readfile(filename):
    with open(filename, 'r+') as f:
        return f.read()


setup(name='tilepy',
      version='1.2.0',
      description='Gravitational waves follow-up scheduling and simulation for IACTs',
      install_requires=[
          'gammapy', #==0.9',
          'astropy', #>=4.1',
          'scipy', #',
          'healpy', #==1.12.9',
          'ipython', #',
          'matplotlib', #==2.2.2',
          'MOCpy', #==0.10.0',
          'numpy', #==1.19.5',
          'pandas',
          'pytz', #==2021.1',
          'ephem', #==3.7.6.0',
          'gdpyc',
          'ligo.skymap~=1.0.7',
      ],
      packages=find_packages(),
      # tests_require=['pytest'],
      author='Monica Seglar-Arroyo ',
      author_email='mosear2@gmail.com ',
      url='https://drf-gitlab.cea.fr/multimessenger-IRFU/cta/gw-follow-up-simulations',
      long_description=readfile('README.md'),
      include_package_data=True,
      package_data={'': ['dataset/*.all']},
      # license='',
      classifiers=[],
      scripts=[],
      )
