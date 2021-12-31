
# 1. Overview


## Purpose
This file is intended to provide an overview of the COTS Star Tracker.

IT IS STRONGLY RECOMMENDED THAT YOU READ THE ENTIRETY OF THE /DOCS DIRECTORY AND THE SCRIPTS WITHIN /EXAMPLES BEFORE USING THE SOFTWARE.


## Table of Contents
* [1.1 Use](#1.1-use)
* [1.2 Background](#1.2-background)
* [1.3 Overview](#1.3-overview)
* [1.4 High-Level Description of Use](#1.4-high-level-description-of-use)
* [1.5 System Testing](#1.4-system-testing)
* [1.6 Further Reading](#1.5-further-reading)
* [1.7 Contact](#1.6-contact)


## 1.1 Use
It's recommended that users start here in order to understand how to successfully install, configure, and operate the system.
The above topics are broken out into their own Markdown-formatted files to allow for easy navigation.  First-time users are recommended
to start with this document and then proceed through the others in the docs directory in numerical order.


## 1.2 Background
(As of 2020) Star trackers are common sensors in spaceflight that use energy gathered by detector arrays combined with star catalogs to determine
the inertial attitude of the boresight of the detector array.  They typically have an accuracy better than 50 arc seconds and come
in a variety of weights and sizes.  Many have optical heads that fit within a 150 mm by 150 mm by 150 mm volume (not including a baffle)
and have a mass of a few kilograms.  Star tracker typically cost hundreds of thousands of dollars and take many months to build.  
These instruments, like many other avionics, continue to get smaller and cheaper.  CubeSat star trackers typically fit with a
50 mm by 50 mm by 100 mm volume and cost tens of thousands of dollars, but still have relatively long lead times.

The CubeSat spacecraft standard has facilitated the rise of low-cost launch opportunities and a growing capability for these
small spacecraft to meaningfully contribute to science and technology development and demonstration.  These spacecraft are
developed by a variety of groups, ranging from government, to industry, to academia.  All of these spacecraft are extremely
volume constrained and some groups developing CubeSats are also very financially constrained.  This drives the need for further
cost and volume reduction for enabling technologies for these spacecraft, including attitude determination.  That has lead to 
this development effort.


## 1.3 Overview
The COTS Star Tracker is a set of algorithms that's intended to enable a user to get high accuracy attitude information from
a picture of a star field (much like the aforementioned star trackers) taken with a COTS camera.  These algorithms are meant
to be generic and efficient so that they return reasonably accurate results on commonly available, low size, weight, and power (SWAP)
COTS computers (e.g. Raspberry Pi 3B+).  The algorithms are written in Python 3 as it's a highly capable, open source, high-level programming language
that has a large user base, is cross-platform, and features relatively descriptive Traceback errors to accelerate development.

The algorithms are open source, the COTS cameras are relatively low cost and lead time, and the COTS computers are relatively
low cost and lead time as well.  This should result in a capability that's easy to acquire and low cost.  The hope is that this
enables small satellite developers to accelerate their schedules, reduce their costs, and improve their capability.

There are many other documented approaches to star identification that can be replicated in software and even some publicly-available
software repositories that function as attitude estimation utilities operating from images of star fields.  The COTS Star Tracker
system differs from these as it's demonstrated with a variety of hardware (both imagers and computers) and is built to be rapidly
deployed on a variety of systems by those that aren't experts in imaging systems or attitude determination.

The accuracy and efficacy of the COTS Star Tracker system depends on many things, including the selected camera, lens, quality of the calibration,
the parameters used to configure the algorithm, mounting errors, etc.  Even when properly calibrated and configured, the accuracy of the COTS Star
Tracker system is still far below that of commercially available systems (like those mentioned above).  The COTS Star Tracker is not
meant to replace commercially available systems-- it's meant to enhance the capability of spacecraft that aren't able to accommodate them.

## 1.4 High-Level Description of Use
You are again STRONGLY encouraged to read all of the documents and scripts in the /examples directory before using this software.  In order
to provide a big-picture view to help the others docs make more sense in context, a high-level description of the steps required to use this
software and implement your own star tracker system follows:

1. Select a camera and lens
    * Not all cameras are created equal.  We do not recommend those that only support basic UVC drivers
2. Interface the camera to your computer
    * This is both hardware and software.  Example interface scripts for several camera manufacturers are included in the /examples directory.
3. Calibrate your camera and lens
    * Set your camera/lens focus and aperture as they would be in flight (you can use stars on a clear night to help you set these well).
    * Collect a series of images with the camera/lens of either checkerboard targets or stars such that you get good coverage across the field of view.
    * Appropriately process the images (using checkerboard cal for the checkerboard and Tetra cal for the star fields).
    * Calibration may be an iterative process-- double check all resulting values in the .json files for sanity (e.g. your focal length isn't 1E9) and aim for an RMS reprojection error < 1 pixel
4. Build a star catalog
    * This requires the camera calibration file created in the previous step as well as user input to provide the limiting apparent magnitude for the stars that will be included in the catalog.
    * Note that this can take a considerable amount of time on low-resource single-board computers and is best done on a "modern" desktop workstation
5. Attempt to process images of stars.
    * This can be either HWIL or an existing image set.  Tune various parameters until you get the balance of accuracy, sky coverage, and solve time that you need for your application.
    * To assess your configuration's accuracy, you can use the included scripts to process your images through astrometry.net for an independent solution.  Included scripts can also compare the resulting products.
6. Demonstrate your system works end-to-end with hardware-in-the-loop
    * If this wasn't done in the previous step, ensure your entire system works, including all the hardware

## 1.5 System Testing
A preliminary version of the COTS Star Tracker algorithm was tested on a variety of image sets, including some taken from a synthetic generation system,
some taken on-orbit, and some taken terrestrially.  These image sets were processed on a Raspberry Pi 3B+, an Odroid XU4Q, and a Dell Precision 7720.
The resulting accuracy across all platforms was fairly consistent with a mean error of about 200 arc seconds.  The solve time varied significantly, with
most solve times on the Dell being around 0.25 s, on the Pi being 2 s, and the Odroid being 1.5 s.

System updates and performance analysis is still in work.


## 1.6 Further Reading
For more information on a preliminary version of the system, please review the following paper: https://ntrs.nasa.gov/citations/20200001376


## 1.7 Contact
For questions/comments/concerns, please do not hesitate to get in touch and we will reply as able.

Please direct all correspondence to the COTS Star Tracker Point of Contact, Sam Pedrotty: samuel.m.pedrotty@nasa.gov



