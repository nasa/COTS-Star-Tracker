#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
install.py

This program is intended to install all
prerequisite packages for the star tracker
software on a computer running Windows

Note that this is useful for data processing
and analysis, but may not work for hardware-in-
the-loop testing.

'''

################################
#LOAD LIBRARIES
################################
import os
import sys
import time

################################
#MAIN CODE
################################
# install/update python stuff
os.system('pip3 install pip --upgrade --user')
os.system('pip3 install opencv-contrib-python --user')
os.system('pip3 install psutil --user')
os.system('pip3 install imageio') #required for catalog creation
os.system('pip3 install astropy') #required for catalog creation
os.system('pip3 install pandas') #required for catalog creation
os.system('pip3 install statistics --user')
os.system('pip3 install astroquery') #required for astrometry verification

#must install freetype2 dev pkg first??
os.system('pip3 --no-cache-dir install matplotlib --user')
os.system('pip3 install --upgrade setuptools --user')
os.system('pip3 --no-cache-dir install scipy --user')

#install module
home = os.getcwd()
os.chdir('..')
os.chdir('..')
os.chdir('py_src/star_tracker')
os.system('pip3 install .')
os.chdir(home)

# install/update stuff for IDS cam
os.system('pip3 install pyueye --user')

print("\n\nInstallation complete.  Please restart the computer!") 
print("NOTE: no camera interfaces were installed during this process.  Other scripts/software may have to be run to install camera software interfaces and drivers.")

