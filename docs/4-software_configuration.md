
# 4. Software Configuration

NOTE: it's assumed that the user has already generated a camera calibration file in step 3.1 prior to starting this procedure.


## Purpose
This file is intended to provide detailed information on the configuration of the software for use with the star tracker algorithm.
For information on the use and verification of the software, see 5-software_use.md


## Table of Contents
* [4.1 Star Catalog Creation](#4.1-star-catalog-creation)
* [4.2 Star Tracker Algorithm Setup](#4.2-star-tracker-algorithm-setup)
* [4.3 Troubleshooting](#4.3-troubleshooting)


## 4.1 Star Catalog Creation
This step uses the camera parameters determined in Step 3.1, the Hipparcos star catalog, and user input to create star catalog specific to
the configuration of this system.

1. open the py_src/tools/star_catalog_creator.py script in a text editor
2. update the b_thresh variable to the appropriate value.  This is the brightest apparent magnitude that will be included in the catalog.
3. point the script to the camera parameter file.  Camera parameters are used to more intelligently create the catalog.
4. point the script to the directory where images are stored to be used to create the darkframe
5. Run the star_catalog_creator.py script
    * This can be fairly time consuming on slower single board computers like the Raspberry Pi 3B+
    * After the catalog creation, the script will also attempt to create a darkframe from images in the data directory.  If there are none, it will fail, but the catalog will still be good to go.


## 4.2 Star Tracker Algorithm Setup
This step describes how to create a script that captures images, processes them, and returns an attitude estimate.  It's strongly recommended that users review
the script(s) in the 'examples' directory before continuing.

1. Create a new Python 3 script
2. Using the camera interface previously demonstrated, call the routines required to configure the camera and capture a picture
3. Pass the following information to the main.star_tracker function in the star_tracker library
    * first argument (string): image filename
    * second argument (string):camera parameter filename
    * m keyword argument (numpy data): m file contents
    * q keyword argument (numpy data): q file contents
    * x_cat keyword argument (numpy data): x_cat file contents
    * k keyword argument (numpy data): k file contents
    * indexed_star_pairs keyword argument (numpy data): indexed star pair file contents
    * Optionally, you can also provide:
        * darkframe_file keyword argument (string): the name of the darkframe file
            * this defaults to 'None', in which case a darkframe will not be used.
        * undistort_img_bool keyword argument (bool): the flag to enable or disable distortion compensation
            * this defaults to 'True'.
        * n_stars keyword argument (integer): the number of stars that will be used to determine catalog matches (always starting with the brightest)
            * this defaults to 30.  Fewer stars may reduce solve times, but may reduce sky coverage.
        * isa_thresh keyword argument (float): inner star angle threshold
            * this defaults to 0.0008, which is equivalent to about 5px.  You may want to provide a smaller number.
        * nmatch keyword argument (integer): minimum number of catalog matches required to estimate attitude
            * this defaults to 6.  It's not recommended to use a lower number as it may result in incorrect attitude estimates.
        * low_thresh_pxl_intensity keyword argument (integer): the lower threshold for pixel intensity
            * this defaults to 'None', which will then automatically determine a lower value.  If centroiding is failing when stars are obviously visible in the image, manually setting this to a lower value may help.
        * hi_thresh_pxl_intensity keyword argument (integer): the upper threshold for pixel intensity
            * this defaults to 'None', which will then automatically be set to the highest pixel value in the image.
        * min_star_area keyword argument (integer): the minimum pixel area to be considered a star (pixels)
            * this defaults to 4 pixels.
        * max_star_area keyword argument (integer): the maximum pixel area to be considered a star (pixels)
            * this defaults to 36 pixels.  If you notice that some large stars are not being considered as centroids, try increasing this number.
        * watchdog (float): the maximum time to be spent in the triangle_isa_id function, essentially bounding solve time.
            * this defaults to None, which applies a timeout duration of 3600 seconds.  Typically, the solve time is bimodal, with one peak being correct solutions and the other being failures to solve.  If you have failures to solve taking up lots of compute time, set this to prevent the long solve times that almost always result in failures to solve.
        * graphics keyword argument (bool): the flag to enable or disable graphical output
            * this defaults to 'False', which disables graphical output.  When troubleshooting, set to 'True' to see visualizations of the algorithm throughout the process.
        * verbose keyword argument (bool): the flag to enable or disable verbose printing
            * this defaults to 'False', which disables verbose printing.  When troubleshooting, set to 'True' to see text output of the algorithm throughout the process.
4. The star_tracker function will return a quaternion estimate, list of matches, number of matches, x_obs, and the image after manipulation
    * first return (array): quaternion
        * quaternion is right-handed, scalar last.
    * second return (array): idmatch
        * array of matches from the triangle_isa_id function
    * third return (integer): nmatches
        * this is the number of successful matches to the catalog of the solution used to estimate the attitude
    * fourth return (array): x_obs
        * vector observations that were the inputs to the triangle_isa_id function
    * fifth return (image array): img
        * the input image after manipulation

Example call:

    q_est, idmatch, nmatches, x_obs, rtrnd_img = main.star_tracker(
            "my_star_image.jpg", "my_cam_file.json", m=m, q=q, x_cat=x_cat, k=k, indexed_star_pairs=indexed_star_pairs, 
            nmatch=7, n_stars=25, graphics=True, verbose=True)


## 4.2 Troubleshooting
This step described generic troubleshooting steps to help overcome errors in configuring the software and getting it to work for the first time.

This software has been designed with the intention to compensate for bad inputs and provide detailed errors when it fails to do so.  That being said, the software
is still very immature and there are many exceptions that have yet to be caught.  This software attempts to provide errors and warnings:

* Warnings: these are situations where the software encounters an issue, but can overcome it by making an assumption.  These assumptions should be understood by the user so they are not left "silent."
    * A warning should be displayed in the terminal output in this way: "WARNING [star_tracker.function_name]: description of the warning"
    * Warnings don't necessarily have to be acted on, but the user should be aware when they are happening and what their impact is
* Errors: these are situations where the software can not continue and will exit.  This requires a user to make a change to the software inputs to allow the software to function.
    * An error should be displayed in the terminal output in this way: "ERROR [star_tracker.function_name]: description of the error"