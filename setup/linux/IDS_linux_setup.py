#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
ids_linux_setup.py

This program is intended to download and
install latest IDS cam drivers for linux

'''

################################
#LOAD LIBRARIES
################################
import os

################################
#MAIN CODE
################################
# install/update stuff for IDS cam
os.system('sudo pip3 install pyueye')
os.system('sudo apt-get --yes install python3-qt4 python3-qt4-doc')

print("To install the IDS drivers, you must:")
print("  * download the right driver from: https://en.ids-imaging.com/download-ueye-emb-hardfloat.html ")
print("       (sadly, you have to create an account)")
print("  * extract the tar file")
print("  * run the appropriate .run file")

