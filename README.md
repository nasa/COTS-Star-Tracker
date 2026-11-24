
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
* [Related NASA Resources](#related-nasa-resources)

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
1. images taken of a known checkerboard calibration target or star field images combined with a 3rd party tool are used to computer camera and lens parameters
2. camera parameters are used (along with user input) to create a reference catalog from the Hipparcos star database
3. based on user-input parameters, images are processed to identify star centroids
4. based on user-input parameters, the star centroids are matched to the reference catalog
5. if a suitable match is found, a quaternion representing the attitude of the camera boresight is returned


## Questions Comments and Concerns

The best way to get these things addressed is to open an issue on GitHub (https://github.com/nasa/COTS-Star-Tracker) and/or email the developer (samuel.m.pedrotty@nasa.gov).
Don't hesitate to reach out-- it's likely others also have your questions/comments/concerns.  We encourage users to add 
capability and fix bugs and then submit merge/pull requests to have their improvements brought back to the repository for all to benefit from.


## High-Level Release Notes

* v1.0: initial release
* v1.1: improved docs, bug fixes, new camera calibration method
* v1.2: improved docs, bug fixes, addition of error codes


## Related NASA Resources

### Additional Information and Flight Information

* Initial publication of COTS Star Tracker System: https://ntrs.nasa.gov/citations/20200001376
* Flown on R5 Spacecraft 2 (R5-S2) and R5 Spacecraft 4 (R5-S4): https://digitalcommons.usu.edu/smallsat/2025/RA-S4-2025/1/
* Flown on R5 Spacecraft 9 (R5-S9) and R5 Spacecraft 10 (R5-S10): https://digitalcommons.usu.edu/smallsat/2026/all2026/9/

### GIANT

The COTS Star Tracker is "flight" software for performing attitude estimation from star field images and geometric camera calibration.  
For other image-based navigation needs we encourage you to check out our partner project, the [Goddard Image Analysis and Navigation Tool (GIANT)](https://github.com/nasa/giant), which is ground software providing a larger set of navigation, calibration, and science capabilities.


### Orion Optical Navigation Image Processing Software (OpNAV)

The Orion emergency return system utilizing optical navigation was demonstrated on Artemis I.
A key component is this software, which processes the images and produces navigation measurements. 
The software is an application within the Core Flight Software of the Orion camera controller module and may be requested here: https://software.nasa.gov/software/MSC-26456-1



