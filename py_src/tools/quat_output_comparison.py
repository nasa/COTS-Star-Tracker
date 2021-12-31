#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
quat_output_comparison.py

This script is intended to load two
different files that contain quaternions
and determine the difference between
corresponding quaternions
'''


################################
#IMPORT MODULES
################################
import os
import time
import pandas as pd
import matplotlib.pyplot as plt
from star_tracker.support_functions import *
from star_tracker.array_transformations import *


################################
#USER INPUT
################################
enable_plots = True #enable/disable plotting
plot_filename_prefix = '' #the prefix for your plot filenames
plot_title_prefix = '' #the prefix for your plot titles
file1_name = '' #the name of the first file to compare
file2_name = '' #the name of the second file to compare

file1_image_name_fieldname = 'image name'
file1_quat_scalar_fieldname = 'qs'
file1_quat_vec1_fieldname = 'qv0'
file1_quat_vec2_fieldname = 'qv1' 
file1_quat_vec3_fieldname = 'qv2'
file1_solvetime_fieldname = 'image solve time (s)'

file2_image_name_fieldname = 'image name'
file2_quat_scalar_fieldname = 'qs'
file2_quat_vec1_fieldname = 'qv0'
file2_quat_vec2_fieldname = 'qv1'
file2_quat_vec3_fieldname = 'qv2'
file2_solvetime_fieldname = 'image solve time (s)'


################################
#MAIN CODE
################################

rad2deg = 180/pi

# load data
print("loading data")
df1 = pd.read_csv(file1_name, sep=',', index_col=False)
df2 = pd.read_csv(file2_name, sep=',', index_col=False)

# get total number of file 1 solns
file1_solns = 0
file1_soln_names = []
file1_solvetime_good = []
file1_solvetime_bad = []
for n in range(0,len(df1[file1_image_name_fieldname])):
    if df1[file1_quat_scalar_fieldname][n] < 999:
        file1_solns+=1
        file1_soln_names+=[df1[file1_image_name_fieldname][n]]
        file1_solvetime_good+=[df1[file1_solvetime_fieldname][n]]
    else:
        file1_solvetime_bad+=[df1[file1_solvetime_fieldname][n]]


# get total number of file 2 solns
file2_solns = 0
file2_soln_names = []
file2_solvetime_good = []
file2_solvetime_bad = []
for n in range(0,len(df2[file2_image_name_fieldname])):
    if df2[file2_quat_scalar_fieldname][n] < 999:
        file2_solns+=1
        file2_soln_names+=[df2[file2_image_name_fieldname][n]]
        file2_solvetime_good+=[df2[file2_solvetime_fieldname][n]]
    else:
        file2_solvetime_bad+=[df2[file2_solvetime_fieldname][n]]


print("\nFile 1 ("+file1_name+") has "+str(file1_solns)+" solutions of "+str(len(df1[file1_quat_scalar_fieldname]))+" total ("+str((file1_solns/len(df1[file1_quat_scalar_fieldname]))*100)[:5]+"%)")
print("\nFile 2 ("+file2_name+") has "+str(file2_solns)+" solutions of "+str(len(df2[file2_quat_scalar_fieldname]))+" total ("+str((file2_solns/len(df2[file2_quat_scalar_fieldname]))*100)[:5]+"%)")


#compare and identify number of images with no matching solution
print("processing data")
common_names = []
theta_err = []
x_err = []
y_err = []
z_err = []
nonzero = 0

for name in file1_soln_names:
    for file2_image_name in file2_soln_names:
            if (file1_image_name == file2_image_name) or (file1_image_name in file2_image_name) or (file2_image_name in file1_image_name):
                common_names+=[file2_image_name]
                q1_s = df1.loc[df1[file1_image_name_fieldname] == name, file1_quat_scalar_fieldname].values[0]
                q1_v0 = df1.loc[df1[file1_image_name_fieldname] == name, file1_quat_vec1_fieldname].values[0]
                q1_v1 = df1.loc[df1[file1_image_name_fieldname] == name, file1_quat_vec2_fieldname].values[0]
                q1_v2 = df1.loc[df1[file1_image_name_fieldname] == name, file1_quat_vec3_fieldname].values[0]
                att_est1 = quat2attitude_matrix([q1_v0,q1_v1,q1_v2], q1_s)

                q2_s = df2.loc[df2[file2_image_name_fieldname] == name, file2_quat_scalar_fieldname].values[0]
                q2_v0 = df2.loc[df2[file2_image_name_fieldname] == name, file2_quat_vec1_fieldname].values[0]
                q2_v1 = df2.loc[df2[file2_image_name_fieldname] == name, file2_quat_vec2_fieldname].values[0]
                q2_v2 = df2.loc[df2[file2_image_name_fieldname] == name, file2_quat_vec3_fieldname].values[0]
                att_est2 = quat2attitude_matrix([q2_v0,q2_v1,q2_v2], q2_s)

                if q1_s == q2_s and q1_v0 == q2_v0 and q1_v1 == q2_v1 and q1_v2 == q2_v2:
                    theta_err += [0]
                    x_err += [0]
                    y_err += [0]
                    z_err += [0]
                else:
                    phi, theta = attitude_error(att_est1, att_est2)
                    theta_err += [theta*rad2deg]
                    x_err += [float(phi[0])*rad2deg]
                    y_err += [float(phi[1])*rad2deg]
                    z_err += [float(phi[2])*rad2deg]
                    nonzero+=1



print("\n     "+str(len(common_names))+" common solutions identified between the files")
print("         .... "+str((nonzero/len(common_names))*100)+"% were not identical")


# TODO save output


# plot
if enable_plots:
    print("Plotting...")

    n=0
    plt.figure(n)
    plt.hist(df1[file1_solvetime_fieldname], bins='auto')
    plt.ylabel('#')
    plt.xlabel('solve time(s)')
    plt.title(plot_title_prefix+' file 1 solve time (s)')
    plt.savefig(plot_filename_prefix+'time1_overall.jpg')

    n+=1
    plt.figure(n)
    plt.hist(file1_solvetime_good, bins='auto')
    plt.ylabel('#')
    plt.xlabel('solve time(s)')
    plt.title(plot_title_prefix+' file 1 successful solve time (s)')
    plt.savefig(plot_filename_prefix+'time1_success.jpg')

    n+=1
    plt.figure(n)
    plt.hist(file1_solvetime_bad, bins='auto')
    plt.ylabel('#')
    plt.xlabel('solve time(s)')
    plt.title(plot_title_prefix+' file 1 unsuccessful solve time (s)')
    plt.savefig(plot_filename_prefix+'time1_fail.jpg')

    n+=1
    plt.figure(n)
    plt.hist(df2[file2_solvetime_fieldname], bins='auto')
    plt.ylabel('#')
    plt.xlabel('solve time(s)')
    plt.title(plot_title_prefix+' file 2 solve time (s)')
    plt.savefig(plot_filename_prefix+'time2_overall.jpg')

    n+=1
    plt.figure(n)
    plt.hist(file2_solvetime_good, bins='auto')
    plt.ylabel('#')
    plt.xlabel('solve time(s)')
    plt.title(plot_title_prefix+' file 2 successful solve time (s)')
    plt.savefig(plot_filename_prefix+'time2_success.jpg')

    n+=1
    plt.figure(n)
    plt.hist(file2_solvetime_bad, bins='auto')
    plt.ylabel('#')
    plt.xlabel('solve time(s)')
    plt.title(plot_title_prefix+' file 2 unsuccessful solve time (s)')
    plt.savefig(plot_filename_prefix+'time2_fail.jpg')

    n+=1
    plt.figure(n)
    plt.hist(x_err, bins='auto')
    plt.ylabel('#')
    plt.xlabel('error (deg)')
    plt.title(plot_title_prefix+' X delta')
    plt.savefig(plot_filename_prefix+'xerr.jpg')

    n+=1
    plt.figure(n)
    plt.hist(y_err, bins='auto')
    plt.ylabel('#')
    plt.xlabel('error (deg)')
    plt.title(plot_title_prefix+' Y delta')
    plt.savefig(plot_filename_prefix+'yerr.jpg')

    n+=1
    plt.figure(n)
    plt.hist(z_err, bins='auto')
    plt.ylabel('#')
    plt.xlabel('error (deg)')
    plt.title(plot_title_prefix+' Z delta')
    plt.savefig(plot_filename_prefix+'zerr.jpg')

    n+=1
    plt.figure(n)
    plt.hist(theta_err, bins='auto')
    plt.ylabel('#')
    plt.xlabel('error (deg)')
    plt.title(plot_title_prefix+' total delta (deg)')
    plt.savefig(plot_filename_prefix+'thetaerr.jpg')


    plt.show()

print('Done')








