
# COTS Star Tracker README

The commercial off-the-shelf (COTS) camera and computer-based Star Tracker (COTS Star Tracker) is intended
to enable reasonably high accuracy attitude estimation using COTS cameras and computers.  This open source Python 3-based
software is intended to be paired with a user-provided COTS camera and computer.  Several common COTS
single board computers (SBCs) and COTS cameras have been demonstrated with this software through its highly-
automated installation, configuration, and operation.  It's hoped that this software will enable improved mission
performance and lower costs.  It's also hoped that the software's capability and performance will continue to improve
with the help of the community.

Please see the docs directory for more information on the configuration and use of the software.

IT IS STRONGLY RECOMMENDED THAT YOU READ THE ENTIRETY OF THE /DOCS DIRECTORY AND THE SCRIPTS WITHIN /EXAMPLES FIRST.


## Table of Contents

* [High-Level Use](#high-level-use)
* [High-Level Functionality](#high-level-functionality)
* [Questions Comments and Concerns](#questions-comments-and-concerns)
* [High-Level Release Notes](#high-level-release-notes)


## High-Level Use

It's hoped that users will be able to quickly and cheaply add high accuracy attitude determination to their spacecraft with the following steps:
1. acquire suitable COTS camera and computer
2. integrate COTS camera and computer
3. move COTS Star Tracker code to computer
4. connect computer to the internet and run the appropriate install script
5. calibrate camera
6. create star catalog
7. verify function with live-sky testing
8. integrate into spacecraft system


## High-Level Functionality

The code functions as outlined below:
1. images taken of a known checkerboard calibration target or starfield images combined with a 3rd party tool are used to computer camera and lens parameters
2. camera parameters are used (along with user input) to create a reference catalog from the Hipparcos star database
3. based on user-input parameters, pictures are processed to identify star centroids
4. based on user-input parameters, the star centroids are matched to the reference catalog
5. if a suitable match is found, a quaternion representing the attitude of the camera boresight is returned


## Questions Comments and Concerns

The best way to get these things addressed is to open an issue on GitHub (https://github.com/nasa/COTS-Star-Tracker) and/or email the developer (samuel.m.pedrotty@nasa.gov).
Don't hesitate to reach out-- it's likely others also have your questions/comments/concerns.  We encourage users to add 
capability and fix bugs and then submit merge/pull requests to have their improvements brought back to the repository for all to benefit from.


## High-Level Release Notes

* v1.0: initial release
