
# Appendix: Function Description

Note: this section is still in work.


## Purpose
This file is intended to provide additional information about core functions of the star tracker algorithm.



**`cam_matrix()`**
##### Input
 - **`dx`**:	dx component of intrinsic camera matrix
 - **`dy`**:	dy component of intrinsic camera matrix
 - **`up`**:	principal u-axis
 - **`vp`**:	principal v-axis
 - **`sk`**:	skew component of the camera matrix
##### Output
 - **`c`**:	intrinsic camera matrix inverse


**`cam_matrix_inv()`**
##### Input
 - **`c`**:	intrinsic camera matrix
##### Output
 - **`cinv`**: inverse camera matrix


**`cam2fov()`**
##### Input
 - **`cinv`**:	inverse camera matrix
 - **`nrow`**:  number of pixel rows in the image
 - **`ncol`**:  number of pixel columns in the image
##### Output
 - **`fov`**: field of view (3x1 array, [diagonal fov, x, y])
 

**`read_cam_json()`**
##### Input
- **`cam_config_file`**: camera configuration data file
- **`cam_resolution`**: resolution parameters of the camera in pixels
- **`pixel_pitch`**: pixel pitch parameters of the camera in millimeters
- **`focal_length`**: focal length of the camera in millimeters
#####Output
- **`camera_matrix`**: intrinsic camera matrix \
$\begin{bmatrix}
   fx & skew & up \\
   0 & fx & up \\
   0 & 0 & 1
\end{bmatrix}$
- **`cam_resolution`**:
- **`dist_coefs`**: distortion coefficients of the camera


**`r[X]_passive()`** 
##### Input
 - **`theta`**:	angle of rotation
##### Output
 - **`r`**: rotation matrix using passive transformation


**`attitude_matrix2quat()`**
##### Input
- **`a`**: rotation matrix from inertial reference frame to that of the s/c camera
##### Output
- **`q`**: quaternion resulting from attitude matrix (first three quantities are vector components, the fourth is the scalar component)


**`quat2attitude_matrix()`**
##### Input
- **`q13`**: vector quantity of quaternion describing rotation from inertial reference frame to s/c camera
- **`q4`**: scalar quantity of quaternion describing rotation from inertial reference frame to s/c camera
##### Output
- **`a`**: rotation matrix from inertial reference frame to that of the s/c camera


**`vector_dot()`**
##### Input
- **`a`**: first vector array (nxm) in dot product operation
- **`b`**: second vector array (nxm) in dot product operation\
(n = number of elements in each vector, m = number of vectors in array)
##### Output
- dot product result (nx1) of the two vector arrays


**`vector_norm()`**
##### Input
- **`v`**: vector array (nxm)\
(n = number of elements in each vector, m = number of vectors in array)
##### Output
- vector norm (nx1) of the vector array

**`normalize_vector_array()`**
##### Input
- **`v`**: vector array (nxm)\
(n = number of elements in each vector, m = number of vectors in array)
##### Output
- normalized vector array (nxm)


**`camera2homogenous()`**
##### Input
- **`p`**: x,y coordinates in camera space
##### Output
- homogeneous coordinates normalized to z = 1 (x/Z, y/Z, 1)


**`camera2vector()`**
##### Input
- **`xy`**: x,y vector array in camera frame (nx2)
##### Output
- array of unit vectors resulting from camera frame to world frame (nx3)


**`pixel2vector()`**
##### Input
- **`cinv`**: inverse camera matrix
- **`p`**: vector array of pixel coordinates (u,v) in image (nx2)
##### Output
- array of unit vectors resulting from pixel space to world frame (nx3)


**`homogeneous2vector()`**
##### Input
- **`cinv`**: inverse camera matrix
- **`h`**: focal-length normalized homogeneous coordinates (x/Z, y/Z, 1)
##### Output
- array of unit vectors resulting from homogeneous space to world frame (nx3)


**`vector2pixel()`**
##### Input
- **`c`**: instrinsic camera matrix
- **`v`**: 3xn vector array of world coordinates
##### Output
2xn vector array of pixel coordinates (u,v)


**`vector2homegeneous()`**
##### Input
- **`c`**: instrinsic camera matrix
- **`v`**: 3xn vector array of world coordinates
##### Output
3xn vector array of camera (x,y) coordinates


**`vector_array_transform()`**
##### Input
- **`transform_matrix`**: transformation matrix (nxn)
- **`v`**: vector array undergoing transformation (nxm)
##### Output
- vector array resulting from matrix multiplication of vector array (nxm)


**`stars_in_fov()`**
##### Input  
 - **`u_i`**: unit vectors of each star in the inertial frame  
 - **`t_i2c`**: transformation matrix from inertial to camera frame  
 - **`c`**: camera parameters matrix  
 - **`nrow`**: number of pixel rows in the image  
 - **`ncol`**: number of pixel columns in the image  
 - **`fov`**: angle describing the circular field of view
##### Output
 - **`star_idx_fov`**: index of stars within the rectangular field of view  


**`basic_cam()`**
##### Input
- **`nrow`**: number of rows
- **`ncol`**: number of columns
- **`dx`**:	dx component of intrinsic camera matrix
- **`dy`**:	dy component of intrinsic camera matrix
- **`up`**:	principal u-axis
- **`vp`**:	principal v-axis
- **`sk`**:	skew component of the camera matrix
##### Output
- **`c`**: intrinsic camera matrix
- **`nrow`**: number of rows
- **`ncol`**: number of columns
- **`fov`**: angle describing the circular field of view


**`read_cam_file()`**
##### Input
- **`cam_config`**: camera configuration data file
- **`data_rows`**: number of rows in the data
- **`data_cols`**: number of columns in the data 
- **`data_path`**: path of the data
##### Output
- **`fields`**: 
- **`cam1filename`**:
- **`cam2filename`**:


 
**`kvector()`**
##### Input  
 - **`inputCat`**: catalog used to create the kvector
##### Output  
 - **`k`**: k-vector containing index values from the star pair interstar angles
 - **`m`**: slope of the interpolation line
 - **`q`**: intercept of the interpolation line
 - **`sorted_cat`**:	catalog of star pair indices matched with interstar angle

**`equatorial2vector()`**
##### Input  
 - **`ra_de`**: right ascension and declination of each star in catalog [nx2 or 2xn]
 - **`axis`**: axis determination \
 (Axis = 0: each column represents individual star RA and DEC\
  Axis = 1: each row represents individual star RA and DEC)
##### Output  
 - **`u`**: line-of-sight unit vectors to the ith star as seen by an observer at the BCRF origin
 
**`vector2equatorial()`**
##### Input
 - **`v`**:
 - **`axis`**: axis determination \
 (Axis = 0: each column represents individual star unit vectors\
  Axis = 1: each row represents individual star unit vectors)
##### Output


**`mas2rad()`**
##### Input
 - **`pm_rade`**: proper motion values in milli-arcseconds
##### Output
- proper motion values in radians

**`lpq_orthonormal_basis()`**
##### Input
 - **`ra_de`**: right ascension and declination of each star in catalog [nx2 or 2xn]
##### Output
- **`l`**: line-of-sight vectors to each star in catalog
- **`p`**: array of cross products of the camera direction with each star in catalog  (z x l, z = transpose([0 0 1])
- **`q`**: array of cross products of line-of-sight vectors and p (l x p)

**`proper_motion_correction()`**
##### Input
- **`ra_de`**:	RA and DE in degrees
- **`pm_rade`**:  proper motion of RA and DE in mas/yr
- **`plx`**: parallax in mas
- **`rB`**:	BCRF position in km of celestial object that the spacecraft orbits
- **`t`**: time of observation
- **`t_ep`**: time of the catalog epoch
##### Output
- **`los`**: array of line-of-sight unit vectors to each star as seen by an observer at the BCRF origin
- **`u`**: array of vectors to each star after correcting for proper motion

**`create_star_catalog()`**
##### Input
 - **`starcat`**: data file containing the Hipparcos catalog  
 - **`brightness_thresh`**: desired magnitude threshold
 - **`col_data`**: 
 - **`col_brightness`**: column in star catalog containing the brightness values
 - **`col_pm`**: catalog column containing the proper motion values
 - **`col_par`**: catalog column containing parallax values
 - **`row_start`**: catalog row where data starts
 - **`save_vals`**: option to save values or just return them from function
 - **`fov`**: angle describing the circular field of view
 - **`rB`**: BCRF position of the celestial body the spacecraft is orbiting
 - **`cat_ep`**: catalog epoch (in years using Julian calendar)
##### Output
 - **`k`**: k-vector containing index values from the star pair interstar angles
 - **`m`**: slope of the interpolation line
 - **`q`**: intercept of the interpolation line
 - **`d_cat`**: catalog of star pair indices matched with interstar angles
 - **`u`**: unit vectors of stars in catalog in ICRF
 - **`d_cat_idx`**: catalog of star pair indices matched with interstar angles including indices of all star pairs in catalog
 
**`read_dir_images()`**
##### Input
- **`folder`**: path to images (may be relative or absolute path)
- **`numImages`**: number of images used in creating the dark frame
##### Output
- **`images`**: nxmxr array containing r images, each of size nxm
- **`img_dtype`**: type of data used in images (uint8, uint16, etc.)

**`create_darkframe()`**
##### Input
- **`folder`**: path to images (may be relative or absolute path)
- **`numImages`**: number of images used in creating the dark frame
##### Output
- **`dark_frame`**: image used to subtract noise in background

 

**`nchoosek()`**
##### Input  
 - **`n`**:  number of options in list
 - **`k`**:  number of items in combination
##### Output  
 - **`r`**: unique combinations of list
 
 **`interstar_angle()`**
##### Input  
 - **`star_pair`**: coordinates of two stars in a pair  
##### Output  
 - **`theta`**: angle between two stars
 
**`enhanced pattern shifting()`**
##### Input  
 - **`candidateIdx`**: ordered candidate star indices from list
##### Output  
 - **`pKernel`**: patterned kernel for output into star ID algorithm  

**`kvec_values()`**
##### Input
- **`inputCat`**: input catalog
- **`s_idx`**: indices
##### Output
- **`m`**: slope of the interpolation line
- **`q`**: intercept of the interpolation line
 
**`ksearch()`**
##### Input  
 - **`kVec`**:  k-vector containing index values from the star pair interstar angles
 - **`xin`**:  search value of k-vector
 - **`d_thresh`**: range within which star pairs are considered a match 
 - **`d_cat`**:	catalog of star pair indices matched with interstar angle
 - **`m`**: slope of the interpolation line
 - **`q`**: intercept of the interpolation line
##### Output   
 - **`matchid`**: indices of all candidates from d_cat falling within specified range

**`fullobsmatch()`**
##### Input
 - **`x_obs`**:  unit vectors from image pixel coordinates
 - **`x_cat`**:  unit vectors of all stars in catalog that fall within constraints
 - **`d_thresh`**: range within which star pairs are considered a match 
##### Output
 - **`idmatch`**:  indices of star pair matches
 - **`nmatches`**: number of star matches from image
 
 **`attitude_svd()`**
##### Input
 - **`ei`**: catalog of unit vectors in the inertial frame, 3xn matrix  
 - **`es`**: measured unit vectors in the sensor frame, 3xn matrix  
##### Output
 - **`T`**: rotation matrix from inertial frame to sensor frame  

**`find_candidate_stars()`**
##### Input
- **`img`**: image containing star field
- **`low_thresh`**: low threshold for finding centroids, used in findContours() and connectedComponentsWithStats()
- **`hi_thresh`**: high threshold for finding centroids, used in findContours() and connectedComponentsWithStats()
##### Output
- **`found_centroids`**: centroids found using findContours()
- **`cand_centroids`**: centroids that have enough area to categorize as candidate stars

**`star_tracker()`**
##### Input
- **`img_file`**: path to image
- **`cam_config_file`**: camera configuration file
- **`darkframe_file`**: dark frame file
- **`nmatch`**: number of matches
- **`x_cat`**: 
- **`d_cat_idx`**:
- **`d_thresh`**: 
- **`k`**: k-vector containing all star pairs in catalog
- **`m`**: slope of the interpolation line
- **`q`**: intercept of the interpolation line
- **`low_thresh_pxl_intensity`**: low threshold for finding centroids, used in findContours() and connectedComponentsWithStats()
- **`hi_thresh_pxl_intensity`**: high threshold for finding centroids, used in findContours() and connectedComponentsWithStats()
##### Output
- **`q_est`**: estimated attitude quaternion found from star tracker algorithm
- **`nmatches`**: number of matches from quaternion in the catalog in the image
- **`x_obs`**: candidate star unit vectors
- **`camera_matrix`**: intrinsic camera matrix


