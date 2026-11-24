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

#pull down and install ximea drivers
os.system('sudo wget  https://updates.ximea.com/public/ximea_linux_sp_beta.tgz')
os.system('sudo tar xzf ximea_linux_sp_beta.tgz')
os.chdir('package')
os.system('sudo ./install')

print("Ximea camera software/driver installation complete!")
print(" If you encounter issues running the software, you may have to add utf8 headers to Python API in order to get them to function:  -*- coding: utf-8 -*- ")
