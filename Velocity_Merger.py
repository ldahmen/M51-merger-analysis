"""
    Adapted:            5/11/2020
    Last update:        5/11/2020
    Author:             Linnea Dahmen, Adapted from Jose Flores Velazquez
    Jose's website:     https://jaf12.github.io/joseflores.github.io/
"""


################################################################################

""" This program illustrates interaction of Galaxies M51 and NGC 5195, allowing
    the user to manipulate initial conditions of the neighboring galaxy's
    velocity. """

################################################################################

import numpy as np
import glob
import sys
import os
import copy
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import patches
from matplotlib.lines import Line2D

#...Set up font size:
fontsize = 14

################################################################################

#---- UNITS OF CODE

"""
    These units make it possible for gravitational constant to equal 1 (G=1).

1 unit of Mass     = 2x10^(10) Solar Masses
1 unit of Distance = 500 pc
1 unit of Time     = 1.2 million years
1 unit of Velocity = 400 km/s
"""

################################################################################

#--- Written by Jose, this function moves the simulation forward.
def acceleration_stars(spiral_coordinates, neighbor_coordinates,
    star_coordinates_x, star_coordinates_y, star_coordinates_z, params, t):

    #---- Unpacking Variables
    X1 = spiral_coordinates[0]
    Y1 = spiral_coordinates[1]
    Z1 = spiral_coordinates[2]

    X2 = neighbor_coordinates[0]
    Y2 = neighbor_coordinates[1]
    Z2 = neighbor_coordinates[2]

    x = star_coordinates_x
    y = star_coordinates_y
    z = star_coordinates_z

    G  = params[0]
    sf = params[1]
    M1 = params[2]
    M2 = params[3]

    #---- Determining Radial Distance (Denominator of Acceleration)
    r1 = np.sqrt((X1-x)**2 + (Y1-y)**2 + (Z1-z)**2 + sf**2)
    r2 = np.sqrt((X2-x)**2 + (Y2-y)**2 + (Z2-z)**2 + sf**2)

    #---- Acceleration for each star
    acceleration_x = (((G*M1)/(r1**3))*(X1-x)) + (((G*M2)/(r2**3))*(X2-x))
    acceleration_y = (((G*M1)/(r1**3))*(Y1-y)) + (((G*M2)/(r2**3))*(Y2-y))
    acceleration_z = (((G*M1)/(r1**3))*(Z1-z)) + (((G*M2)/(r2**3))*(Z2-z))


    return acceleration_x, acceleration_y, acceleration_z, r1, r2



#--- Written by Jose, this function moves the simulation forward.
def acceleration_nuclei(spiral_coordinates, neighbor_coordinates, params, t):

    #---- Unpacking Variables
    X1 = spiral_coordinates[0]
    Y1 = spiral_coordinates[1]
    Z1 = spiral_coordinates[2]

    X2 = neighbor_coordinates[0]
    Y2 = neighbor_coordinates[1]
    Z2 = neighbor_coordinates[2]

    G  = params[0]
    sf = params[1]
    M1 = params[2]
    M2 = params[3]

    #---- Determining Radial Distance (Denominator of Acceleration)
    s = np.sqrt((X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2 + sf**2)

    #---- Acceleration for each nuclei
    Acceleration_1x = ((G*M2)/(s**3)) * (X2-X1)
    Acceleration_1y = ((G*M2)/(s**3)) * (Y2-Y1)
    Acceleration_1z = ((G*M2)/(s**3)) * (Z2-Z1)

    Acceleration_2x = ((G*M1)/(s**3)) * (X1-X2)
    Acceleration_2y = ((G*M1)/(s**3)) * (Y1-Y2)
    Acceleration_2z = ((G*M1)/(s**3)) * (Z1-Z2)

    Acceleration_Spiral      = np.array([Acceleration_1x, Acceleration_1y, Acceleration_1z])
    Acceleration_Neighboring = np.array([Acceleration_2x, Acceleration_2y, Acceleration_2z])

    return Acceleration_Spiral, Acceleration_Neighboring, s


#--- This function is adapted to streamline running the simulation for multiple
#--- Velocities. It takes in a list of velocities, a list of "tags" for saving
#--- the data in a unique way, the number of time steps you'd like to simluate
#--- for, and the number of stars you'd like to populate the simulation.
#--- It uses the functions that Jose has created.
def mult_merg(velocities,vtags,ts,N):
    for p in range(len(velocities)):

        v = velocities[p] # see lines 158 and 159 for use of this variable
        velocity_tag = vtags[p] # this defines the naming scheme of files
        timeStep = ts # see line 175 for the use of this variable
        N = N # see line 191 for the use of this variable

        #---- Galaxy Parameters
        initial_diameter_spiral  = 52.0
        initial_radius_spiral    = initial_diameter_spiral / 2.0
        Number_rings_stars       = 10

        #---- Mass of galaxies
        Spiral_Galaxy_Mass      = 5.0
        Neighboring_Galaxy_Mass = (1.0/4.0) * Spiral_Galaxy_Mass

        #---- Initial Positions of Galaxies
        Spiral_Galaxy_x0 = 0.0
        Spiral_Galaxy_y0 = 0.0
        Spiral_Galaxy_z0 = 0.0

        #Neighboring_Galaxy_x0 = -30.0
        #Neighboring_Galaxy_y0 = -30.0
        Neighboring_Galaxy_x0 = -35.0
        Neighboring_Galaxy_y0 = -50.0
        Neighboring_Galaxy_z0 = 0.0
        Pos_Neigh = np.sqrt(Neighboring_Galaxy_x0**2 + Neighboring_Galaxy_y0**2 + Neighboring_Galaxy_z0**2)

        #---- Initial Velocities of Galaxies
        Spiral_Galaxy_v_x0 = 0.0
        Spiral_Galaxy_v_y0 = 0.0
        Spiral_Galaxy_v_z0 = 0.0

        Neighboring_Galaxy_v_x0 = 0.0
        Neighboring_Galaxy_v_y0 = v
        Neighboring_Galaxy_v_z0 = v

        #---- Softening Factor, Graviational Constant, Degree of Freedom
        #sf  = 3.5
        #sf  = 2.01
        sf = 2.016
        #sf = 1.0
        G   = 1.0
        dof = 3

        #---- Time Step
        dt        = 1.0
        #time_step = 540
        time_step = timeStep
        time      = np.arange(0,time_step+1,dt)
        half_dt    = 0.5 * dt
        half_dt_sq = 0.5 * (dt*dt)

        #---- Parameters
        params = np.array([G, sf, Spiral_Galaxy_Mass, Neighboring_Galaxy_Mass])

        #---- Initial Positions of Stars in Spiral Galaxy

        #...A ring is located every 2.6 distance appart
        radius_each_ring = np.arange(0, initial_radius_spiral+(initial_radius_spiral/Number_rings_stars), initial_radius_spiral/Number_rings_stars)
        radius_each_ring = radius_each_ring[1::]

        #...Degrees of each star in every quadrant
        #degree_each_star_Q1 = np.arange(0, 90.0, 180.0/25.0)
        num = N/500*25
        degree_each_star_Q1 = np.arange(0, 90.0, 180.0/num)
        degree_each_star_Q4 = (degree_each_star_Q1[1::] * -1.04)
        degree_each_star_Q2 = degree_each_star_Q1 + 90
        degree_each_star_Q3 = degree_each_star_Q4 + 270
        degree_each_star_Q4 = degree_each_star_Q4 + 360

        #...Connecting degrees of stars to form circle
        semi_circle_top    = np.concatenate((degree_each_star_Q1,degree_each_star_Q2))
        semi_circle_bottom = np.concatenate((degree_each_star_Q3, degree_each_star_Q4))
        degree_circle      = sorted(np.concatenate((semi_circle_top, semi_circle_bottom)))

        #...Emtpy data containers to store positions of stars in each ring
        x_pos = []
        y_pos = []

        #...Updating empty containers with x and y positions of stars
        for j in np.arange(len(radius_each_ring)):
            for i in np.arange(len(degree_circle)):
                 x_pos.append(radius_each_ring[j]*np.cos(degree_circle[i]*np.pi/180.0))
                 y_pos.append(radius_each_ring[j]*np.sin(degree_circle[i]*np.pi/180.0))

        #...Total Star
        Number_stars                  = len(x_pos)
        Initial_Star_Coordinates      = np.zeros([Number_stars, dof])
        Initial_Star_Coordinates[:,0] = x_pos
        Initial_Star_Coordinates[:,1] = y_pos

        #---- Star Data
        Star_Coordinates_x = np.zeros([time_step,Number_stars])
        Star_Coordinates_y = np.zeros([time_step,Number_stars])
        Star_Coordinates_z = np.zeros([time_step,Number_stars])

        Star_Velocities_x = np.zeros([time_step,Number_stars])
        Star_Velocities_y = np.zeros([time_step,Number_stars])
        Star_Velocities_z = np.zeros([time_step,Number_stars])

        Star_Accelerations_x = np.zeros([time_step,Number_stars])
        Star_Accelerations_y = np.zeros([time_step,Number_stars])
        Star_Accelerations_z = np.zeros([time_step,Number_stars])

        #---- Initializing Star Data
        Star_Coordinates_x[0] = Initial_Star_Coordinates[:,0]
        Star_Coordinates_y[0] = Initial_Star_Coordinates[:,1]
        Star_Coordinates_z[0] = Initial_Star_Coordinates[:,2]

        #---- Establishing Data Arrays
        Spiral_Galaxy_Coordinates = np.zeros([time_step,dof])
        Neighb_Galaxy_Coordinates = np.zeros([time_step,dof])

        Spiral_Galaxy_Velocities  = np.zeros([time_step,dof])
        Neighb_Galaxy_Velocities  = np.zeros([time_step,dof])

        Spiral_Galaxy_Accelerations = np.zeros([time_step,dof])
        Neighb_Galaxy_Accelerations = np.zeros([time_step,dof])

        #---- Initializing Data Arrays
        Spiral_Galaxy_Coordinates[0]  = Spiral_Galaxy_x0, Spiral_Galaxy_y0, Spiral_Galaxy_z0
        Spiral_Galaxy_Velocities[0]   = Spiral_Galaxy_v_x0, Spiral_Galaxy_v_y0, Spiral_Galaxy_v_z0

        Neighb_Galaxy_Coordinates[0]  = Neighboring_Galaxy_x0, Neighboring_Galaxy_y0, Neighboring_Galaxy_z0
        Neighb_Galaxy_Velocities[0]   = Neighboring_Galaxy_v_x0, Neighboring_Galaxy_v_y0, Neighboring_Galaxy_v_z0


        #...Initial Acceleration of Nuclei
        spiral_acc, neighb_acc, s = acceleration_nuclei(Spiral_Galaxy_Coordinates[0], Neighb_Galaxy_Coordinates[0], params,0)


        #---- Star Acceleration
        for i in np.arange(len(Star_Coordinates_z)):
            Star_Accelerations_x[i], Star_Accelerations_y[i], Star_Accelerations_z[i], r1, r2 =  acceleration_stars(Spiral_Galaxy_Coordinates[0], Neighb_Galaxy_Coordinates[0], Star_Coordinates_x[i], Star_Coordinates_y[i], Star_Coordinates_z[i], params, 0)

        #---- Velocity Verlet Method For Calculating Positions, Speeds, & Accelerations

        #...Calculations for Center Buldges Positions, Speeds, & Accelerations
        for t in range(1,time_step):

            #---- Updating Coordinates
            Spiral_Galaxy_Coordinates[t] = Spiral_Galaxy_Coordinates[t-1] + Spiral_Galaxy_Velocities[t-1]*dt + half_dt_sq * Spiral_Galaxy_Accelerations[t-1]
            Neighb_Galaxy_Coordinates[t] = Neighb_Galaxy_Coordinates[t-1] + Neighb_Galaxy_Velocities[t-1]*dt + half_dt_sq * Neighb_Galaxy_Accelerations[t-1]

            #---- Current Acceleration
            acc1, acc2, s = acceleration_nuclei(Spiral_Galaxy_Coordinates[t], Neighb_Galaxy_Coordinates[t], params, time[t])
        #    print "s", acc1, acc2
            #---- Updating Velocities
            Spiral_Galaxy_Velocities[t] = Spiral_Galaxy_Velocities[t-1] + half_dt * (acc1 + Spiral_Galaxy_Accelerations[t-1])
            Neighb_Galaxy_Velocities[t] = Neighb_Galaxy_Velocities[t-1] + half_dt * (acc2 + Neighb_Galaxy_Accelerations[t-1])

            Spiral_Galaxy_Accelerations[t] = acc1
            Neighb_Galaxy_Accelerations[t] = acc2

        #...Calculations for Star's Positions, Speeds, & Accelerations
        for t in range(1,time_step):

            for j in range(Number_stars):

                #---- Updating Coordinates
                Star_Coordinates_x[t][j] = Star_Coordinates_x[t-1][j] + Star_Velocities_x[t-1][j]*dt + half_dt_sq * Star_Accelerations_x[t-1][j]
                Star_Coordinates_y[t][j] = Star_Coordinates_y[t-1][j] + Star_Velocities_y[t-1][j]*dt + half_dt_sq * Star_Accelerations_y[t-1][j]
                Star_Coordinates_z[t][j] = Star_Coordinates_z[t-1][j] + Star_Velocities_z[t-1][j]*dt + half_dt_sq * Star_Accelerations_z[t-1][j]

                #---- Current Acceleration
                acc_x, acc_y, acc_z, r1, r2 = acceleration_stars(Spiral_Galaxy_Coordinates[t], Neighb_Galaxy_Coordinates[t], Star_Coordinates_x[t][j], Star_Coordinates_y[t][j], Star_Coordinates_z[t][j], params, time[t])
        #        print "r's", acc_x, acc_y, acc_z

                #---- Updating Velocities
                Star_Velocities_x[t][j] = Star_Velocities_x[t-1][j] + half_dt * (acc_x + Star_Accelerations_x[t-1][j])
                Star_Velocities_y[t][j] = Star_Velocities_y[t-1][j] + half_dt * (acc_y + Star_Accelerations_y[t-1][j])
                Star_Velocities_z[t][j] = Star_Velocities_z[t-1][j] + half_dt * (acc_z + Star_Accelerations_z[t-1][j])

                Star_Accelerations_x[t][j] = acc_x
                Star_Accelerations_y[t][j] = acc_y
                Star_Accelerations_z[t][j] = acc_z

        #---------------------------------------------------------------------------
        #--- This section of the function creates the radius rings at 20% mass and
        #--- 80%.
        t = time_step
        distances = []

        for k in range(t):
            temp = []

            s_xc = Star_Coordinates_x[k]
            s_yc = Star_Coordinates_y[k]
            s_zc = Star_Coordinates_z[k]

            g_xc = Spiral_Galaxy_Coordinates[k][0]
            g_yc = Spiral_Galaxy_Coordinates[k][1]
            g_zc = Spiral_Galaxy_Coordinates[k][2]

            for i in range(s_xc.size):
                x = s_xc[i]
                y = s_yc[i]
                z = s_zc[i]

                r = np.sqrt((x-g_xc)**2+(y-g_yc)**2 +(z-g_zc)**2)
                temp.append(r)

            temp = np.array(temp)
            temp = np.sort(temp)
            distances.append(temp)
        distances =  np.array(distances)

        r80 = []
        r20 = []

        for i in range(t):
            r80.append(distances[i][int(0.8*s_xc.size)])
            r20.append(distances[i][int(0.2*s_xc.size)])

        r80 = np.array(r80)
        r20 = np.array(r20)


        #---------------------------------------------------------------------------
        #--- Function that generates frames to create video, written by Jose and
        #--- adapted to have plot titles and radius indicating rings.
        def plot_2D(field = 'video', save_frame=False, eps_tag=False):

            #---- Set up save tag:
            plot_tag =  field
            save_tag = plot_tag

            #---- Creating multiple snapshot images
            for i in np.arange(time_step):
                print("Snapshot:", i)

                #--- Plotting
                fig, ax = plt.subplots()
                stars = ax.scatter(Star_Coordinates_x[i],Star_Coordinates_y[i], c= 'c', label = '$Stars$')
                neighb = ax.scatter(Neighb_Galaxy_Coordinates[i][0],Neighb_Galaxy_Coordinates[i][1], label='$NGC$ $5195$',   c='r', s=250)
                spir = ax.scatter(Spiral_Galaxy_Coordinates[i][0],Spiral_Galaxy_Coordinates[i][1], label='$M51$ $Center$', c='k', s=100 )
                circle1 = plt.Circle((Spiral_Galaxy_Coordinates[i][0],Spiral_Galaxy_Coordinates[i][1]),radius=r80[i], color='green',fill=False,label='$80\%$ $Radius$')
                circle2 = plt.Circle((Spiral_Galaxy_Coordinates[i][0],Spiral_Galaxy_Coordinates[i][1]),radius=r20[i], color='purple',fill=False,label='$20\%$ $Radius$')
                ax.add_artist(circle1)
                ax.add_artist(circle2)

                #--- radius indicators for the legend
                green_circle = Line2D([0], [0], marker='o', color='green',
                                fillstyle='none', markersize=15)
                purple_circle = Line2D([0], [0], marker='o', color='purple',
                                fillstyle='none', markersize=10)

                #--- Plotting Labels
                ax.set_xlabel('$X-Pos$ [$500pc$]', fontsize=17)
                ax.set_ylabel('$Y-Pos$ [$500pc$]', fontsize=17)
                ax.set_title('$Merger:$ $N=$ '+ str(N) + ' $and$ $velocity=$ ' + str(v))

                ax.set_xlim(-60,60) # delete these if you like your axes to move
                ax.set_ylim(-50, 50) # ^^       ^^      ^^      ^^      ^^

                #--- create legend with specified labels
                ax.legend((stars, neighb, spir, green_circle, purple_circle), ('$Stars$', '$NGC$ $5195$', '$M51$ $Center$', '$80\%$ $Radius$','$20\%$ $Radius$')
                    ,loc='upper left', scatterpoints = 1)

                #-------------------------------------------------------------------
                #--- This section of the function creates a unique directory and
                #--- stores all of the created images within said directory

                #---- Makes Frames True/False
                if save_frame == False:

                    #...If saved plots directory doesn't exit, make one.
                    if not os.path.exists('./merger2D_'+velocity_tag): os.mkdir('./merger2D_'+velocity_tag)

                    #...Set up file type:
                    file_type = '.pdf'
                    if eps_tag == True: file_type = '.eps'

                    #...Save pdf:
                    plt.savefig('./merger2D_'+velocity_tag+'/'+save_tag+'_'+file_type)

                elif save_frame == True:

                    #...If saved frames directory doesn't exit, make one.
                    if not os.path.exists('./merger2D_'+velocity_tag): os.mkdir('./merger2D_'+velocity_tag)

                    #...Set up file type:
                    file_type = '.png'

                    #...Save pngs:
                    plt.savefig('./merger2D_'+velocity_tag+'/frame'+str(i)+'.png',dpi=300)

                #...Report and close
                print("SAVE: ", save_tag)
                plt.close()

        #---------------------------------------------------------------------------
        #--- Runs the previous function
        #--- Make the frames

        plot_2D(field = 'video', save_frame=True, eps_tag=False)

        #--- Save the positional data of the stars and the galaxy so that we can do
        #--- further analysis
        np.save('./Star_Coordinates_x_'+velocity_tag+'.npy',Star_Coordinates_x)
        np.save('./Star_Coordinates_y_'+velocity_tag+'.npy',Star_Coordinates_y)
        np.save('./Star_Coordinates_z_'+velocity_tag+'.npy',Star_Coordinates_z)
        np.save('./Spiral_Galaxy_Coordinates_'+velocity_tag+'.npy',Spiral_Galaxy_Coordinates)


################################################################################

N = 1000 # was done with 500 originally, 1000 provides better resolution
velocities = [0.2, 0.34, 0.5, 0.75, 1, 1.5, 2, 2.5, 3] #change this to your liking
ts = 540 # to see the initial perturbation, we use ts = 300. Default is 540.

#... create velocity tags
#... current format, i.e. v = 0.2: '0.2_N1000'
vtags = []
for i in range(len(velocities)):
    tag = str(velocities[i]) + '_N' + str(N) # add whatever you want here
    vtags.append(tag)

#... Note velocities and vtags must be the same size
#--- RUN THE FUNCTION
mult_merg(velocities,vtags,ts,N)
