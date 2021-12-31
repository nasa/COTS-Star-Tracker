#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
astrometry_output_comparison.py

This script is intended to load a star tracker
output file and an astrometry.net output file
and then compare their outputs in order to enable
the use of astrometry.net as an independent way to
assess the accuracy of the results.
'''


################################
#IMPORT MODULES
################################
import os
import csv
import time
import numpy as np
import pandas as pd
import output_converter
import matplotlib.pyplot as plt
from math import *
from datetime import datetime
from star_tracker.support_functions import *
from star_tracker.array_transformations import *


################################
#USER INPUT
################################
enable_plots = True #enable/disable plotting
plot_filename_prefix = '' #the prefix for your plot filenames
plot_title_prefix = '' #the prefix for your plot titles

star_tracker_filename = '' #the name of the first file to compare, output from the star tracker
star_tracker_image_name_fieldname = 'image name'
star_tracker_quat_scalar_fieldname = 'qs'
star_tracker_quat_vec1_fieldname = 'qv0'
star_tracker_quat_vec2_fieldname = 'qv1' 
star_tracker_quat_vec3_fieldname = 'qv2'
star_tracker_solvetime_fieldname = 'image solve time (s)'

astrometry_filename = '' #the name of the second file to compare, output from astrometry
astrometry_image_name_fieldname = 'image_name'
astrometry_right_asc_fieldname = 'ra'
astrometry_dec_fieldname = 'dec'


################################
#MAIN CODE
################################

rad2deg = 180/pi
deg2rad = pi/180

# load data
print("loading data")
st_df = pd.read_csv(star_tracker_filename)
astro_df = pd.read_csv(astrometry_filename)

# get total number of file 1 solns
file1_solns = 0
file1_soln_names = []
file1_solvetime_good = []
file1_solvetime_bad = []
for n in range(0,len(st_df[star_tracker_image_name_fieldname])):
    if st_df[star_tracker_quat_scalar_fieldname][n] < 999:
        file1_solns+=1
        file1_soln_names+=[st_df[star_tracker_image_name_fieldname][n]]
        file1_solvetime_good+=[st_df[star_tracker_solvetime_fieldname][n]]
    else:
        file1_solvetime_bad+=[st_df[star_tracker_solvetime_fieldname][n]]


# get total number of file 2 solns
file2_solns = 0
file2_soln_names = []
file2_solvetime_good = []
file2_solvetime_bad = []
for n in range(0,len(astro_df[astrometry_image_name_fieldname])):
    if astro_df[astrometry_right_asc_fieldname][n] < 999:
        file2_solns+=1
        file2_soln_names+=[astro_df[astrometry_image_name_fieldname][n]]



print("\nFile 1 ("+star_tracker_filename+") has "+str(file1_solns)+" solutions of "+str(len(st_df[star_tracker_quat_scalar_fieldname]))+" total ("+str((file1_solns/len(st_df[star_tracker_quat_scalar_fieldname]))*100)[:5]+"%)")
print("\nFile 2 ("+astrometry_filename+") has "+str(file2_solns)+" solutions of "+str(len(astro_df[astrometry_right_asc_fieldname]))+" total ("+str((file2_solns/len(astro_df[astrometry_right_asc_fieldname]))*100)[:5]+"%)")


#compare and identify number of images with no matching solution
print("\nprocessing data...")
common_names = []
theta_err = []
st_ra = []
st_dec = []
astro_ra = []
astro_dec = []
q1_s = []
q1_v0 = []
q1_v1 = []
q1_v2 = []
nonzero = 0
r=1

for st_name in file1_soln_names:
    the_name = st_name.split('/')
    st_comp_name = the_name[-1].split('\\')[-1]

    for astro_name in file2_soln_names:

            the_name = astro_name.split('/')
            astro_comp_name = the_name[-1].split('\\')[-1]

            if st_comp_name == astro_comp_name:

                common_names+=[st_comp_name]

                # extract ST quat
                q1_s += [st_df.loc[st_df[star_tracker_image_name_fieldname] == st_name, star_tracker_quat_scalar_fieldname].values[0]]
                q1_v0 += [st_df.loc[st_df[star_tracker_image_name_fieldname] == st_name, star_tracker_quat_vec1_fieldname].values[0]]
                q1_v1 += [st_df.loc[st_df[star_tracker_image_name_fieldname] == st_name, star_tracker_quat_vec2_fieldname].values[0]]
                q1_v2 += [st_df.loc[st_df[star_tracker_image_name_fieldname] == st_name, star_tracker_quat_vec3_fieldname].values[0]]

                # convert ST quat to RA/DEC
                euler_matrix = np.zeros([1,3])
                euler_matrix[0] = output_converter.conversion.convert_quaternion('ZXZ', q1_v0[-1], q1_v1[-1], q1_v2[-1], q1_s[-1], degrees=True)[2]
                euler_matrix[0,0] = euler_matrix[0,0] - 90
                euler_matrix[0,1] = 90 - euler_matrix[0,1]
                euler_matrix[0,2] = euler_matrix[0,2] + 180

                #print("------------------------")
                #print(euler_matrix[:,0])
                #print(euler_matrix[:,1])
                #print(euler_matrix[:,2])

                st_ra+= list(euler_matrix[:,0])
                st_dec+= list(euler_matrix[:,1])
                #print(st_ra)
                #print(st_dec)

                # convert ST RA/Dec to unit vector
                x_st = r*cos(st_dec[-1]*deg2rad)*cos(st_ra[-1]*deg2rad)
                y_st = r*cos(st_dec[-1]*deg2rad)*sin(st_ra[-1]*deg2rad)
                z_st = r*sin(st_dec[-1]*deg2rad)

                #extract astrometry RA/Dec
                astro_ra += [astro_df.loc[astro_df[astrometry_image_name_fieldname] == astro_name, astrometry_right_asc_fieldname].values[0]]
                astro_dec += [astro_df.loc[astro_df[astrometry_image_name_fieldname] == astro_name, astrometry_dec_fieldname].values[0]]
                #wrap vals to more closely align with the star tracker's output
                if astro_ra[-1] > 180: astro_ra[-1]=(astro_ra[-1]-360)
                if astro_dec[-1] > 180: astro_dec[-1]=(astro_ra[-1]-360)
                #print(astro_ra)
                #print(astro_dec)

                # convert astrometry RA/Dec to unit vector
                x_astro = r*cos(astro_dec[-1]*deg2rad)*cos(astro_ra[-1]*deg2rad)
                y_astro = r*cos(astro_dec[-1]*deg2rad)*sin(astro_ra[-1]*deg2rad)
                z_astro = r*sin(astro_dec[-1]*deg2rad)

                # calculate dot product
                a_dot_b = x_st*x_astro+y_st*y_astro+z_st*z_astro
                a_mag = sqrt(x_st**2+y_st**2+z_st**2)
                b_mag = sqrt(x_astro**2+y_astro**2+z_astro**2)

                #calculate delta angle
                theta_err += [acos(a_dot_b/(a_mag*b_mag))*rad2deg]


print("\n     "+str(len(common_names))+" common solutions identified between the files\n")
print("\n...processing complete!")

# save output
the_data = {'image_name':common_names,"ST_qs":q1_s,"ST_qv0":q1_v0,"ST_qv1":q1_v1,"ST_qv2":q1_v2,"ST_RA_deg":st_ra,"ST_Dev_deg":st_dec,"Astro_RA_deg":astro_ra,"Astro_Dec_deg":astro_dec,"delta_angle_deg":theta_err}
now = str(datetime.now())
now = now.split('.')
now = now[0]
now = now.replace(' ','_')
now = now.replace(':','-')

#write data
keys=sorted(the_data.keys())
with open(os.path.join(os.getcwd(),  now+'_astrometry_compare.csv'),'w', newline='') as csv_file:
    writer=csv.writer(csv_file)
    writer.writerow(keys)
    writer.writerows(zip(*[the_data[key] for  key in keys]))

# plot
if enable_plots:
    print("Plotting...")

    n=0
    plt.figure(n)
    plt.hist(st_df[star_tracker_solvetime_fieldname], bins='auto')
    plt.ylabel('#')
    plt.xlabel('solve time(s)')
    plt.title(plot_title_prefix+' star tracker solve time (s)')
    plt.savefig(now+'_'+plot_filename_prefix+'time1_overall.jpg')

    n+=1
    plt.figure(n)
    plt.hist(file1_solvetime_good, bins='auto')
    plt.ylabel('#')
    plt.xlabel('solve time(s)')
    plt.title(plot_title_prefix+' star tracker successful solve time (s)')
    plt.savefig(now+'_'+plot_filename_prefix+'time1_success.jpg')

    n+=1
    plt.figure(n)
    plt.hist(file1_solvetime_bad, bins='auto')
    plt.ylabel('#')
    plt.xlabel('solve time(s)')
    plt.title(plot_title_prefix+' star tracker unsuccessful solve time (s)')
    plt.savefig(now+'_'+plot_filename_prefix+'time1_fail.jpg')

    n+=1
    plt.figure(n)
    plt.hist(np.array(theta_err), bins='auto')
    plt.ylabel('#')
    plt.xlabel('error (deg)')
    plt.title(plot_title_prefix+' total delta (deg)')
    plt.savefig(now+'_'+plot_filename_prefix+'thetaerr.jpg')

    n+=1
    plt.figure(n)
    plt.plot(np.array(st_ra),'o', label = "RA")
    plt.plot(np.array(st_dec),'o', label = "Dec")
    plt.ylabel('angle (deg)')
    plt.xlabel('image')
    plt.title(plot_title_prefix+' star tracker right ascension and declination (deg)')
    plt.legend(loc="best")
    plt.savefig(now+'_'+plot_filename_prefix+'st_ra_dec.jpg')

    n+=1
    plt.figure(n)
    plt.plot(np.array(q1_s),'o')
    plt.plot(np.array(q1_v0),'o')
    plt.plot(np.array(q1_v1),'o')
    plt.plot(np.array(q1_v2),'o')
    plt.ylabel('quat components')
    plt.xlabel('image')
    plt.title(plot_title_prefix+' star tracker quat')
    plt.savefig(now+'_'+plot_filename_prefix+'st_quat.jpg')

    n+=1
    plt.figure(n)
    plt.plot(np.array(astro_ra),'o', label = "RA")
    plt.plot(np.array(astro_dec),'o', label = "Dec")
    plt.ylabel('angle (deg)')
    plt.xlabel('image')
    plt.title(plot_title_prefix+' astrometry right ascension and declination (deg)')
    plt.legend(loc="best")
    plt.savefig(now+'_'+plot_filename_prefix+'astro_ra_dec.jpg')

    plt.show()

print('Done')








