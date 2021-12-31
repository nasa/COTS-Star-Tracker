
# 2. Software Install


## Purpose
This file is intended to provide detailed instructions on how to install the COTS Star Tracker software.


## Table of Contents
* [2.1 Camera Interface](#2.1-camera-interface)
* [2.2 COTS Star Tracker Install](#2.2-cots-star-tracker-install)


## 2.1 Camera Interface
This step implements the requisite command and control software for the user-selected camera.

It's assumed that user will select a camera and computer to integrate into a star tracker system.
It's also assumed that the camera will have a Python3 API.  While this isn't required, it allows
for the most seamless integration with the star tracker algorithm.  It's up to the user to install
all software/drivers/etc. required for the operation of their camera.  Some installation and 
interface code is included for Ximea, IDS, and Raspberry Pi cameras in the /examples directory.


## 2.2 COTS Star Tracker Install
This step installs the star tracker software onto the target computer.

1. Copy the COTS Star Tracker code to the target computer
    * This can either be a direct clone or a standard copy operation
2. Connect the target computer to the internet
    * The install scripts assume the system is connected to the internet, although the file updates can be done offline with additional effort on the user's part
3. Run the appropriate installer script
    * For Ubuntu operating systems, use the setup/linux/ubuntu/install.py script.  Be sure to run it as sudo:  sudo python3 install.py
    * For Raspbian operating systems, use the setup/linux/raspbian/install.py script.  Be sure to run it as sudo: sudo python3 install.py
    * For Windows operating systems (which are intended to be used as troubleshooting/development only), use the setup/windows/install.py script
4. Restart the target computer

At this point, all dependencies should be installed and the star tracker code should be able to be called as a Python 3 module.

