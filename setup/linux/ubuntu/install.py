#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
install.py

This program is intended to install all
prerequisite packages for the star tracker
software on a computer running Ubuntu

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

# Ensure locale and time have been set
usr_in = input("\nBefore running this script, the system locale and time should be set.  Has this been done? [Y/N]: ")
usr_in=str(usr_in)
usr_in=usr_in.lower()
if usr_in != 'y' and usr_in != 'yes' and usr_in != 'heck yeah':
    print("\n\nPlease set the system locale and time, restart the system, and then try again.")
    sys.exit()

# Ensure system is up to date
usr_in = input("\nBefore running this script, the system/package managers should be updated to ensure it has access to the latest/best package repositories.  Has this been done? [Y/N]: ")
usr_in=str(usr_in)
usr_in=usr_in.lower()
if usr_in != 'y' and usr_in != 'yes' and usr_in != 'heck yeah':
    print("\n\nPlease update the system/package managers, and then try again.")
    sys.exit()

# Install/update system packages
print("\n\n Installing system packages (this requires an internet connection)...\n")
os.system('sudo apt-get --yes install gfortran')
os.system('sudo apt-get --yes install libblas-dev') # scipy req
os.system('sudo apt-get --yes install python3-pip')

# Define paths
ubuntu_install_dir = os.getcwd()
linux_dir = os.path.dirname(ubuntu_install_dir)
setup_dir = os.path.dirname(linux_dir)
requirements_full_filename = os.path.join(setup_dir,"requirements.txt")
top_dir = os.path.dirname(setup_dir)
star_tracker_src_dir = os.path.join(top_dir, "py_src", "star_tracker")

# Install Python packages
print("\n\n Installing Python packages (this requires an internet connection)...\n")
os.system("pip3 --no-cache-dir install -r "+requirements_full_filename) # --no-cache-dir helps avoid low memory errors
print("\n\n ...complete.\n\n")

# Install star tracker package
print("\n\n Installing Star Tracker package...\n")
os.chdir(star_tracker_src_dir)
os.system('pip3 install .')
os.chdir(ubuntu_install_dir)
print("\n\n ...complete.\n\n")


print("Installation complete (note: no camera interfaces were installed during this process).  The computer may need to be restarted.\n\n")

