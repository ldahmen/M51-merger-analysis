"""
    Created:            5/11/2020
    Last update:        5/11/2020
    Author:             Linnea Dahmen, to be used with an Velocity_Merger.py and
                        Isolated_M51.py.
    Jose's website:     https://jaf12.github.io/joseflores.github.io/
"""


################################################################################

""" This program takes positional data from the merger simulation and does
    an analysis based on the perturbation of M51. """

################################################################################
import numpy as np
import glob
import sys
import os
import copy
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
import matplotlib.font_manager
from numpy import trapz

#--- This function calculates the ratio of interest: the ratio of r80/r20 radii
#--- compared between the isolated galaxy and the merger.
#--- it takes in the velocity tag, corresponding to the directory containing
#--- the data, and returns an array of a ratio calculation at every time step.
def ratio_calc(vtag):
    velocity_tag = vtag

    Star_Coordinates_x = np.load('../sim_data_300s/Star_Coordinates_x_'+velocity_tag+'.npy')
    Star_Coordinates_y = np.load('../sim_data_300s/Star_Coordinates_y_'+velocity_tag+'.npy')
    Star_Coordinates_z = np.load('../sim_data_300s/Star_Coordinates_z_'+velocity_tag+'.npy')
    Spiral_Galaxy_Coordinates = np.load('../sim_data_300s/Spiral_Galaxy_Coordinates_'+velocity_tag+'.npy')
    Star_Coordinates_x_ISO = np.load('../sim_data_300s/Star_Coordinates_x_ISO.npy')
    Star_Coordinates_y_ISO = np.load('../sim_data_300s/Star_Coordinates_y_ISO.npy')
    Star_Coordinates_z_ISO = np.load('../sim_data_300s/Star_Coordinates_z_ISO.npy')
    Spiral_Galaxy_Coordinates_ISO = np.load('../sim_data_300s/Spiral_Galaxy_Coordinates_ISO.npy')

    #-------------------------------------------------------------------------------
    ## Create the data for plotting r80/r100 over time
    # Read in all the positional data for every star and the spiral galaxy at a given time stamp.
    # Next, calulate the distance of each star from the spiral galaxy.
    # Order these distances numerically.
    # Place them into a larger array that will hold all distance values for all stars at all time stamps.
    # Finally, for each time stamp, select the 100th (20 percent) and 400th (80 percent) star for the ratio.

    distances = []
    distances_ISO = []
    #using the first 100 time steps
    for k in range(300):
        temp = []
        temp_ISO = []

        s_xc = Star_Coordinates_x[k]
        s_yc = Star_Coordinates_y[k]
        s_zc = Star_Coordinates_z[k]

        g_xc = Spiral_Galaxy_Coordinates[k][0]
        g_yc = Spiral_Galaxy_Coordinates[k][1]
        g_zc = Spiral_Galaxy_Coordinates[k][2]

        s_xc_ISO = Star_Coordinates_x_ISO[k]
        s_yc_ISO = Star_Coordinates_y_ISO[k]
        s_zc_ISO = Star_Coordinates_z_ISO[k]

        g_xc_ISO = Spiral_Galaxy_Coordinates_ISO[k][0]
        g_yc_ISO = Spiral_Galaxy_Coordinates_ISO[k][1]
        g_zc_ISO = Spiral_Galaxy_Coordinates_ISO[k][2]

        for i in range(s_xc.size):
            x = s_xc[i]
            y = s_yc[i]
            z = s_zc[i]

            r = np.sqrt((x-g_xc)**2+(y-g_yc)**2 +(z-g_zc)**2)
            temp.append(r)

        for j in range(s_xc_ISO.size):
            x = s_xc_ISO[i]
            y = s_yc_ISO[i]
            z = s_zc_ISO[i]

            r = np.sqrt((x-g_xc_ISO)**2+(y-g_yc_ISO)**2 +(z-g_zc_ISO)**2)
            temp_ISO.append(r)


        temp = np.array(temp)
        temp = np.sort(temp)
        distances.append(temp)

        temp_ISO = np.array(temp_ISO)
        temp_ISO = np.sort(temp_ISO)
        distances_ISO.append(temp_ISO)
    distances =  np.array(distances)
    distances_ISO = np.array(distances_ISO)

    ratios = []
    '''
    for i in range(300):
        r_80 = distances[i][int(0.8*s_xc.size)]
        r_20 = distances[i][int(0.2*s_xc.size)]
        value = r_80/r_20
        ratios.append(value)
    '''
    for l in range(300):
        r_80 = distances[l][int(0.8*s_xc_ISO.size)]
        r_20 = distances[l][int(0.2*s_xc_ISO.size)]
        value = r_80/r_20

        r_80_ISO = distances_ISO[l][int(0.8*s_xc_ISO.size)]
        r_20_ISO = distances_ISO[l][int(0.2*s_xc_ISO.size)]
        value_ISO = r_80_ISO/r_20_ISO

        ratios.append(value/value_ISO)

    #data for the plot
    ratios=np.array(ratios)
    return(ratios)

#--- This creates frames to turn into a video to compare with the merger.
#--- There is a time ticker that indicates the time of merger.
def plot_frames(velocities,velocity_tags,N):
    for j in range(len(velocities)):
        v = velocities[j]
        velocity_tag = velocity_tags[j]
        print('VELOCITY: ' +velocity_tag)
        ratios = ratio_calc(velocity_tag)
        area = trapz(ratios, dx=1)
        textstr = ' $AUC:$ ' + str('%.3f'%area)
        x = np.arange(len(ratios))
        time_step = len(ratios)

        for i in range(time_step):
            print("Snapshot:", i)

            plt.clf()
            plt.rcParams.update({'font.size': 14})
            plt.gcf().subplots_adjust(bottom=0.15)
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
            plt.text(200, 13, textstr, fontsize=14, verticalalignment='top',bbox =dict(boxstyle="square", fc="w"))
            plt.plot(i,0.5,color='c', marker='v',ms=10)
            plt.plot(x,ratios,color='r',alpha = 0.5)
            plt.scatter(x, ratios, s=5, color='c',alpha = 1)
            plt.xlabel('time ($1.2\:Myr$ $s^{-1}$)')
            plt.ylabel('$[r80/r20]_{merger}/[r80/r20]_{isolated}$')
            plt.ylim(0,14.2)
            plt.title('Perturbation ratio for velocity of ' +str(v) + ', N=' + str(N))

            if not os.path.exists('../ratios_'+velocity_tag): os.mkdir('../ratios_'+velocity_tag)

            #...Set up file type:
            file_type = '.png'

            #...Save pngs:
            plt.savefig('../ratios_'+velocity_tag+'/frame'+str(i)+'.png',dpi=300)

            plt.close()

def single_plot(velocities,velocity_tags,N):
    for i in range(len(velocities)):

        v = velocities[i]
        velocity_tag = velocity_tags[i]
        print('VELOCITY: ' +velocity_tag)
        ratios = ratio_calc(velocity_tag)

        area = trapz(ratios, dx=1)
        textstr = ' $AUC:$ ' + str('%.3f'%area)
        x = np.arange(len(ratios))

        plt.clf()
        plt.rcParams.update({'font.size': 14})
        plt.gcf().subplots_adjust(bottom=0.15)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        plt.text(200, 13, textstr, fontsize=14, verticalalignment='top',bbox =dict(boxstyle="square", fc="w"))
        plt.plot(x,ratios,color='r',alpha = 0.5)
        plt.scatter(x, ratios, s=5, color='c',alpha = 1)
        plt.xlabel('time ($1.2\:Myr$ $s^{-1}$)')
        plt.ylabel('$[r80/r20]_{merger}/[r80/r20]_{isolated}$')
        plt.ylim(0,14.2)
        plt.title('Perturbation ratio for velocity of ' +str(v) + ', N=' + str(N))

        #...Save pngs:
        plt.savefig('ratio_'+velocity_tag+'.png',dpi=300)
        plt.close()
#-------------------------------------------------------------------------------
velocities = [0.2, 0.34, 0.5, 0.75, 1, 1.5, 2, 2.5, 3]
vtags = ['FIX_0.2_N1000', 'FIX_0.34_N1000', 'FIX_0.5_N1000', 'FIX_0.75_N1000', 'FIX_1_N1000',  'FIX_1.5_N1000','FIX_2_N1000', 'FIX_2.5_N1000', 'FIX_3_N1000' ]
N = 1000
plot_frames(velocities,vtags,N)
#single_plot(velocities,vtags,N)
