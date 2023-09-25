#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
setup.py

This script is intended to enable installation

'''

################################
#LOAD LIBRARIES
################################
from setuptools import setup, find_packages

################################
#MAIN CODE
################################
setup(name='COTS_star_tracker',
      version='1.1',
      description='Open-source framework intended to allow inexperienced users to use COTS cameras and computers to estimate the inertial attitude of the captured images.',
      url='https://github.com/nasa/COTS-Star-Tracker',
      author='NASA',
      author_email='samuel.m.pedrotty@nasa.gov',
      license='BSD 3-Clause',
      packages=find_packages(),
      zip_safe=False)
