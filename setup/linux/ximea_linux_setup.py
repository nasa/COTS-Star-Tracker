#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
ximea_linux_setup.py

This program is intended to download and
install latest ximea cam drivers for linux

'''

################################
#LOAD LIBRARIES
################################
import os

################################
#MAIN CODE
################################
#update USB buffer size and install additional bits
os.system('sudo apt-get update && sudo apt-get install build-essential linux-headers-"$(uname -r)"')
os.system('sudo tee /sys/module/usbcore/parameters/usbfs_memory_mb >/dev/null <<<0')

#required for their streamer apps (not necessarily python API?)
#os.system('sudo apt-get install libgtk2.0-dev')
#os.system('sudo apt-get install libgstreamer1.0-dev libgstreamer-plugins-base1.0-dev')

#pull down and install ximea drivers
os.system('sudo wget https://www.ximea.com/downloads/recent/XIMEA_Linux_SP.tgz')
os.system('sudo tar xzf XIMEA_Linux_SP.tgz')
os.chdir('package')
os.system('sudo ./install')

print("Ximea camera software/driver installation complete!")
print(" If you encounter issues running the software, you may have to add utf8 headers to Python API in order to get them to function:  -*- coding: utf-8 -*- ")
