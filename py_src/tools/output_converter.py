#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
Script to import image names and quaternion calculations from provided file
and calculate respective euler angle/spherical coordinates for each image.
It then compiles all of the data into a dictionary which is keyed with the
image name and valued with the quaternion and euler angle estimations.

    *loads data file
    *Imports image names and quaternion estimations
    *Calculates euler angles from quaternions
    *Compiles data into a dictionary
        **Dictionary key: Image file name
        **Dictionary values: Quaternion and Euler Angle Estimations

'''

################################
#LOAD LIBRARIES
################################
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation as R



################################
#USER INPUT
################################
filename = 'your_filename_here' #include full path
image_name_col = 'image name' #the name of the column containing image names
q_scalar_col = 'qs' #the name of the column containing the scalar component of the quaternion
q_vec1_col = 'qv0' #the name of the column containing the first part of the vector component of the quat
q_vec2_col = 'qv1' #the name of the column containing the second part of the vector component of the quat
q_vec3_col = 'qv2' #the name of the column containing the third part of the vector component of the quat



################################
#SUPPORT FUNCTIONS
################################
class conversion:
    def convert_euler(sequence,x,y,z,degrees):
        '''
        Parameters
        ----------
        sequence : String value that dictates the valid rotation order of the
        Euler Angle input (e.g. 'zxz', 'xyz')
        x : Angle of rotation about the x-axis
        y : Angle of rotation about the y-axis
        z : Angle of rotation about the z-axis
        degrees : Boolean value used to indicate whether or not the angle
        inputs were in degrees(degrees=True) or radians (degrees=False)
        
        Returns
        -------
        r_e : Input values in Euler Angle form per the specified rotation sequence
        r_q : Input values in quaternion form s.t. q = xi + yj + zk + w and
        the "quaternion" output is given by [x,y,z,w]
        r_m : Input values in 3x3 rotation matrix form
        '''

        r = R.from_euler(sequence, [y,x,z], degrees)
        r_q = r.as_quat()
        r_m = r.as_matrix()
        r_e = r.as_euler(sequence,degrees)
        return r_q, r_m, r_e
        
    def convert_quaternion(sequence, x, y, z, w, degrees):
        '''
        Parameters
        ----------
        sequence : String value that dictates the valid rotation order of the
        Euler Angle input (e.g. 'zxz', 'xyz')
        --- Quaternion vector does NOT need to be of the unit vector form ---
        x : Coefficient 'x' from the quaternion form: xi + yj + zk + w
        y : Coefficient 'y' from the quaternion form: xi + yj + zk + w
        z : Coefficient 'z' from the quaternion form: xi + yj + zk + w
        w : Coefficient 'x' from the quaternion form: xi + yj + zk + w
        degrees :   Boolean value used to indicate whether or not the Euler
        angle output will be in degrees (degrees=True) or radians (degrees=False)
        
        Returns
        -------
        r_e : Input values in Euler Angle form per the specified rotation sequence
        r_q : Input values in quaternion form s.t. q = xi + yj + zk + w and
        the "quaternion" output is given by [x,y,z,w]
        r_m : Input values in 3x3 rotation matrix form
        '''

        r = R.from_quat([x,y,z,w])
        r_q = r.as_quat()
        r_m = r.as_matrix()
        r_e = r.as_euler(sequence,degrees)
        return r_q, r_m, r_e
    
    
    def convert_matrix(matrix):
        '''
        Parameters
        ----------
        matrix : 3x3 rotation matrix

        Returns
        -------
        r_q : Quaternion equivalent to the input rotation matrix
        '''

        r = R.from_matrix(matrix)
        r_q = r.as_quat()
        return r_q
            
        

################################
#MAIN CODE
################################

if __name__ == "__main__":
    #Load data
    print("\n\nLoading file...")
    df_in = pd.read_csv(filename)
    #print(df_in)
    print("...file loaded!\n\n")

    image_names = df_in[image_name_col]


    ###############################################################################
    # Storing starting quaternion coeeficients in separate matrices

    print("Populating structures...")
    the_length = len(df_in[q_scalar_col])
    q_matrix_x = df_in[q_vec1_col].to_numpy()
    q_matrix_y = df_in[q_vec2_col].to_numpy()
    q_matrix_z = df_in[q_vec3_col].to_numpy()
    q_matrix_w = df_in[q_scalar_col].to_numpy()

    euler_matrix = np.zeros([the_length,3])
    print("...done!")

    ###############################################################################
    # Converting initial quaternions to matrices, transforming their coordinates, outputing them as euler angle sequence

    print("Math...")
    for i in range(the_length):
    
        if (q_matrix_x[i] == 999) or (q_matrix_y[i] == 999) or (q_matrix_z[i] == 999) or (q_matrix_w[i] == 999):
            euler_matrix[i,:] = np.array([0,0,0])

        else:
            euler_matrix[i] = conversion.convert_quaternion('ZXZ', q_matrix_x[i], q_matrix_y[i], q_matrix_z[i], q_matrix_w[i], degrees=True)[2]
            euler_matrix[i,0] = euler_matrix[i,0] - 90
            euler_matrix[i,1] = 90 - euler_matrix[i,1]
            euler_matrix[i,2] = euler_matrix[i,2] + 180
    print("...math done!")

    ###############################################################################
    # Stores data in a pandas data frame and saves it in an excel file

    df = pd.DataFrame((image_names,q_matrix_x,q_matrix_y,q_matrix_z ,q_matrix_w ,euler_matrix[:,0],euler_matrix[:,1],euler_matrix[:,2]),('File','Q,x','Q,y','Q,z','Q,w','3','1','3'),(np.linspace(1,len(image_names),len(image_names)))).T
    #print(df)
    df.to_csv(filename+"_converted_output.csv")

    print("conversion complete!")



