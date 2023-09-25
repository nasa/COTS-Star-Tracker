
# 5. Software Use


## Purpose
This file is intended to provide detailed instructions on implementing and tuning the performance of the star tracker algorithm.


## Table of Contents
* [5.1 Testing and Verifying](#5.1-testing-and-verifying)
* [5.2 Flight Setup](#5.2-flight-setup)
* [5.3 Considerations](#5.2-considerations)


## 5.1 Testing and Verifying
This step describes common methods and built-in tools to help assess the integrated COTS Star Tracker system performance.  It should
be noted that it's up to the user to determine the appropriate level of rigor to verify the performance of their system.

1. Collect imagery from a live sky with a flight-like camera hardware-in-the-loop (HWIL) setup
    * It's a best practice to then analyze these real-time, but it can be done later as well
    * Note that there are some example camera interface files in the /examples directory that can serve as a starting point
2. Save calculated attitudes from the algorithm in a format that allows the attitudes to be matched up with their corresponding image
    * Again, note the output quaternion is right-handed, scalar last.
3. Utilize a verified tool to process the gathered images and create a "truth" reference to which the system performance can be performed
    * The COTS Star Tracker can be used with astrometry.net to create a reference "truth" data set.
        * create an account on astrometry.net and create an API key.
        * place the key in tools/astrometry_process.py
        * add the path to the directory with the images you'd like to process into the script
        * run the script to automatically have astrometry.net solve the images and record the results: sudo python3 astrometry_process.py
        * an interruption in your internet connection may cause the script to fail.  If that happens:
            * delete the results files with the erroneous "failures"
            * set the 'crash_recovery' variable to True
            * restart the astrometry_process.py script
        * note that it may take on average 20 seconds to successfully solve an image and up to 10 minutes to return a failed solution.  The script has a default timeout of 3600 seconds.
        * the script will create a csv file with the name CURRENT_DATE_astrometry_output_FINAL.csv that contains the Right Ascension, Declination, solve time, and image name
            * for those images that failed to solve, the Right Ascension and Declination will be 999
4. Compare the results
    * In order to compare to astrometry.net, use the tools/astrometry_output_comparison.py script
        * edit the USER INPUT section in tools/astrometry_output_comparison.py
        * run the script: python3 astrometry_output_comparison.py
        * when finished, plots should appear comparing the output of the COTS Star Tracker algorithm and astrometry.net
5.  Adjust the camera and/or star tracker parameters and repeat
    * The graphics and verbose flags of the star tracker algorithm can provide significant insight into how the algorithm is processing imagery.  Use the output to refine the function inputs until the desired results are achieved.


## 5.2 Flight Setup
This step describes a possible way to integrate the star tracker algorithm into a flight software system.  It should
be noted that it's up to the user to determine their sensor performance and flight software integration requirements.

1. Take an image when commanded by a master function or take an image at a set rate
2. Immediately process the image (and optionally save it as well)
3. Return the resulting attitude to the appropriate function or process
4. Repeat


## 5.3 Considerations
This section describes some common considerations that are intended to provoke thought about the configuration of the flight system.

1. If you're planning to save images, what are you going to do with them?
    * Do you even have enough bandwidth to downlink them?
2. What are other performance parameters that would be useful to downlink (e.g. CPU usage, RAM usage, thermal information, etc.)?
3. Does the run rate of the COTS Star Tracker algorithm need to be limited to manage computer temperature?


