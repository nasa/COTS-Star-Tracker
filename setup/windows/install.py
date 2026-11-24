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

# Define paths
win_install_dir = os.getcwd()
setup_dir = os.path.dirname(win_install_dir)
requirements_full_filename = os.path.join(setup_dir,"requirements.txt")
top_dir = os.path.dirname(setup_dir)
star_tracker_src_dir = os.path.join(top_dir, "py_src", "star_tracker")

# Install Python packages
print("\n\n Installing Python packages (this requires an internet connection)...\n")
os.system("pip --no-cache-dir install -r "+requirements_full_filename)
print("\n\n ...complete.\n\n")

# Install star tracker package
print("\n\n Installing Star Tracker package...\n")
os.chdir(star_tracker_src_dir)
os.system('pip3 install .')
os.chdir(win_install_dir)
print("\n\n ...complete.\n\n")


print("Installation complete (note: no camera interfaces were installed during this process).  The computer may need to be restarted.\n\n")

