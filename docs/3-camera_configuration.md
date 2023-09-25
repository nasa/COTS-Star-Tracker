
# 3. Camera Configuration

NOTE: it's assumed that the user has already implemented and demonstrated the camera interface in step 2.1 prior to starting this procedure.


## Purpose
This file is intended to provide detailed information on camera configuration and calibration for use with the star tracker algorithm.


## Table of Contents
* [3.1 Camera Calibration](#3.1-camera-calibration)
* [3.2 Camera Configuration](#3.2-camera-configuration)


## 3.1 Camera Calibration
This step involves estimating various camera parameters.  This may be done by collecting a specific set of imagery with the camera and processing them.

The camera calibration step is intended to capture the set of camera parameters required for compensating for distortion and processing the images.   The COTS Star
Tracker uses the Brown camera model.  Note that camera calibration files created by this software assume a principle point fixed to the center of the detector.
This information is stored as JSON files with the below format:
{
  "camera_matrix": [
    [
      1901.7890625,
      0.0,
      908.7976001542993
    ], 
    [
      0.0,
      1900.397216796875,
      594.3062648404884
    ],
    [
      0.0,
      0.0,
      1.0
    ]
  ],
  "dist_coefs": [
    [
      -0.381661276964259,
      0.4990315571567245,
      0.004867363480130749,
      0.006705504855027461,
      -0.9249236176696716
    ]
  ],
  "resolution": [
    1936,
    1216
  ],
  "CameraModel": "Brown",
  "k1": -0.381661276964259,
  "k2": 0.4990315571567245,
  "k3": -0.9249236176696716,
  "p1": 0.004867363480130749,
  "p2": 0.006705504855027461,
  "fx": 1901.7890625,
  "fy": 1900.397216796875,
  "up": 908.7976001542993,
  "vp": 594.3062648404884,
  "skew": 0.0
}

### 3.1.1.a ChArUco Calibration (recommended)
The charuco_cam_cal.py script in the py_src/tools/camera_calibration/checkerboard directory generates an appropriately formatted JSON file when provided with images of a ChArUco target.
Using the charuco_cam_cal.py script is not required, but an appropriately formatted camera configuration JSON file is required.  At this time, users are encouraged to use the
charuco_cam_cal.py script for camera calibration.

To use the run_tetra_cal.py script:
1. Capture images of stars with your camera set in its final configuration (focus, aperture, etc.)
    * Ensure you can discern more than 5 stars in each photo and capture enough photos to ensure there are stars across the entire field of view
1. Update the 'USER INPUT' section of run_tetra_cal.py
    * Point the script to the location of the aforementioned images, optionally including a darkframe.
    * If you know your camera's field of view, update that as well.  Otherwise, leave the default value
2. Capture a series of images (at least 10) where you can discern at least 5 stars in each picture.
    * Move the camera in between each capture so the resulting image set contains a variety of constellation orientations.
    * Ensure that there are no areas in the field of view that don't contain stars in at least one image.
    * This calibration method may fail on images where many stars are clearly visible.  For that reason, be sure to take 3-4x more images than what you need (e.g. 3-4x more than 10) to ensure you have at least 10 images that Tetra has solved for the calibration.
3. Run run_tetra_cal.py: python3 run_tetra_cal.py
    * Note that this does not have to be done on the target computer and can be done "offline" on a faster one.
4. After completing, the script will either output that it has failed or that is has succeeded and some information about the calibration result, including the RMS error.
5. If the script fails, fails to solve more than 50% of the input images or if the RMS error is more than 1 pixel, take another set of adjust the settings (like FoV) in the run_tetra_cal.py file.
    * If you previously captured many more images than the minimum you need, then you may remove the images that failed to solve from the set and re-solve (assuming you have at least 10 that did solve).
6. If the user does not indicate the cal was a failure, the script will then prompt the user for a name for the camera cal file
    * If no name is selected, it will default to generic_cam_params.json
    * The camera cal file will be located in data/cam_config within the star_tracker module directory


### 3.1.1.b Tetra-based Calibration
The run_tetra_cal.py script in the py_src/tools/camera_calibration/tetra directory generates an appropriately formatted JSON file when provided with images of a starfield.
Using the run_tetra_cal.py script is not required, but an appropriately formatted camera configuration JSON file is required.  If a ChArUco target is unavailable, this calibration
method may be the next best option.

To use the run_tetra_cal.py script:
1. Capture images of stars with your camera set in its final configuration (focus, aperture, etc.)
    * Ensure you can discern more than 5 stars in each photo and capture enough photos to ensure there are stars across the entire field of view
1. Update the 'USER INPUT' section of run_tetra_cal.py
    * Point the script to the location of the aforementioned images, optionally including a darkframe.
    * If you know your camera's field of view, update that as well.  Otherwise, leave the default value
2. Capture a series of images (at least 10) where you can discern at least 5 stars in each picture.
    * Move the camera in between each capture so the resulting image set contains a variety of constellation orientations.
    * Ensure that there are no areas in the field of view that don't contain stars in at least one image.
    * This calibration method may fail on images where many stars are clearly visible.  For that reason, be sure to take 3-4x more images than what you need (e.g. 3-4x more than 10) to ensure you have at least 10 images that Tetra has solved for the calibration.
3. Run run_tetra_cal.py: python3 run_tetra_cal.py
    * Note that this does not have to be done on the target computer and can be done "offline" on a faster one.
4. After completing, the script will either output that it has failed or that is has succeeded and some information about the calibration result, including the RMS error.
5. If the script fails, fails to solve more than 50% of the input images or if the RMS error is more than 1 pixel, take another set of adjust the settings (like FoV) in the run_tetra_cal.py file.
    * If you previously captured many more images than the minimum you need, then you may remove the images that failed to solve from the set and re-solve (assuming you have at least 10 that did solve).
6. If the user does not indicate the cal was a failure, the script will then prompt the user for a name for the camera cal file
    * If no name is selected, it will default to generic_cam_params.json
    * The camera cal file will be located in data/cam_config within the star_tracker module directory


### 3.1.1.c Checkerboard Calibration
The checkerboard_cam_cal.py script in the py_src/tools/camera_calibration/checkerboard directory generates an appropriately formatted JSON file when provided with images of a checkerboard target.
Using the checkerboard_cam_cal.py script is not required, but an appropriately formatted camera configuration JSON file is required.  At this time, the results of a checkerboard-based calibration
can vary significantly, so it's only recommended for users familiar with the method that otherwise are unable to use the above methods.

To use the checkerboard_cam_cal.py script:
1. Download and print a suitable checkerboard pattern
    * the OpenCV one works well: https://docs.opencv.org/master/pattern.png
    * be sure your target is flat (without wrinkles or warps)-- a piece of paper by itself is likely insufficient, but it may suffice if taped to a rigid surface like a desk or clip board.
2. Capture a series of images (at least 10) where the entire target is visible in each picture.  Move the target in between each capture so the set contains a variety of orientations and distances from the camera.  Also ensure that there are no areas in the field of view that don't contain a part of the target in at least one image.
3. Move the images to /py_src/tools/camera_calibration/checkerboard and run checkerboard_cam_cal.py: python3 checkerboard_cam_cal.py
    * Note that this does not have to be done on the target computer and can be done "offline" on a faster one.
    * By default, the checkerboard_cam_cal.py script will attempt to process all images in the same directory as it with a .jpg extension
    * If you use a different target than the OpenCV one linked above, you may need to change the number of rows and columns in the script
    * Running this script can be very time consuming on certain lower power computer systems (like a Raspberry Pi 3B+)
4. Once complete, the script will undistort an image in the set and save two versions: cropped (calibresult_cropped.png) and uncropped (calibresult.png) in the same directory as the script.
    * The script will also print some information about the calibration result, including the RMS error.
5. The script will then prompt the user to verify the cal was successful
    * a visual spot check of calibresult_cropped.png and calibresult.png and verifying that the RMS error is less than 1 pixel is considered sufficient to judge if the cal was a success.  The original image file name will also be printed in the terminal.  If the checkerboard is more square, it's a success.  If it's just as bowed as the original or worse, it's a failure.
6. If the user does not indicate the cal was a failure, the script will then prompt the user for a name for the camera cal file
    * If no name is selected, it will default to generic_cam_params.json
    * The camera cal file will be located in data/cam_config within the star_tracker module directory


## 3.2 Camera Configuration
This step involves the appropriate configuration of the camera's parameters.

This is camera dependent and it's recommended to consult the documentation for the camera selected by the user.  The author's experience 
with various cameras has shown that cameras without sufficient driver support (like those with generic UVC drivers) or appropriate configuration
will not result in sufficient imagery of star fields (at least, when taken from the surface of the Earth with favorable weather conditions).

It's strongly encouraged for users to combine knowledge of their spacecraft ConOps and a clear night within minimal light pollution to appropriately
set the gain and exposure levels for their camera.  Note that longer exposure times will limit the attitude rates at which the algorithm can successfully
identify stars and estimate attitude.  It's forward work to determine these limits in certain configurations.





