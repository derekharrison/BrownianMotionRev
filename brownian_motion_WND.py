import numpy as np
import matplotlib.pyplot as plt
from numpy import *
import sys
import matplotlib.animation as manimation
import time

#Elastic collisions with Periodic boundary conditions

#Simulation parameters
N = 150                            #Number of particles in system
R_a = 0.3                          #Radius of particle a
R_b = 0.1                          #Radii of particle b
L = 3.0                            #Size of periodic domain
max_vel = 2.0                      #Max velocity components of the particles
max_t = 1.0                        #Max simulation time
dt = 1e-2                          #Time step size
m_a = 1.0                          #Particle mass a
m_b = 0.4                          #Particle mass b

start_sim = 0
num_sim = 1

err = 1e-6                         #Small number
h = L - R_b - err                  #Domain limit particle injection

#Initialization
E_wall = L
W_wall = -L
N_wall = L
S_wall = -L

E_lab = 0
W_lab = 1
N_lab = 2
S_lab = 3

NE_lab = 4
NW_lab = 5
SW_lab = 6
SE_lab = 7

C_lab = 8

r_wall_x = N + 1
l_wall_x = N + 2
u_wall_y = N + 3
l_wall_y = N + 4

x = np.zeros(N)
y = np.zeros(N)

Ri = np.zeros(N)
Ri[0] = R_a
for i in range(1, N):
    Ri[i] = R_b

mi = np.zeros(N)
mi[0] = m_a
for i in range(1, N):
    mi[i] = m_b

X = np.zeros((N,9))
Y = np.zeros((N,9))
    
V_X = np.zeros((N,9))
V_Y = np.zeros((N,9))

#Start Multisims
for sim in range(start_sim, num_sim):
    start_time = time.time()
    overlap_counter = 0
    
    #Generate initial velocities
    v_x = max_vel*np.random.randn(N)
    v_y = max_vel*np.random.randn(N)

    #Generating initial positions of particles
    x[0] = 0
    y[0] = 0
    #Generating initial positions of particles
    for n in range(1, N):
        print(n)
        x[n] = 2 * h * np.random.uniform(0, 1) - h
        y[n] = 2 * h * np.random.uniform(0, 1) - h
        overlap = True
        while overlap == True:
            overlap = False
            for i in range(0, n):
                dx = x[i] - x[n]
                dy = y[i] - y[n]
                if dx*dx + dy*dy < (Ri[i]+Ri[n])*(Ri[i]+Ri[n]):
                    overlap = True
            if overlap == True:
                x[n] = 2 * h * np.random.uniform(0, 1) - h
                y[n] = 2 * h * np.random.uniform(0, 1) - h


    #Particle tracking info
    x_track = x[0]
    x_track_zero = x[0]
    y_track = y[0]
    y_track_zero = y[0]
    t_track = []
    r2_track = []

    #Start simulation code
    time_t = 0
    frame_counter = 0
    while time_t < max_t:       
        print(time_t)
        #Exporting data to file
        if frame_counter == floor(time_t / dt):
            #track particle zero
            r2 = ((x_track - x_track_zero)*(x_track - x_track_zero) + (y_track - y_track_zero)*(y_track - y_track_zero))
            t_track.append(time_t)
            r2_track.append(r2)
                
            my_file_x = 'BMW_norm_dist_x_' + str(frame_counter) + '.txt'
            file_x = open(my_file_x, "w")
            my_file_y = 'BMW_norm_dist_y_' + str(frame_counter) + '.txt'
            file_y = open(my_file_y, "w")

            for index in range(0, size(x)):
                data_x = str(x[index])
                data_y = str(y[index])            
                file_x.write(data_x)
                file_x.write('\n')
                file_y.write(data_y)
                file_y.write('\n')

            file_x.close()
            file_y.close()
            frame_counter = frame_counter + 1

        #Initialise collision related info
        coll_partner_1 = 0
        coll_partner_2 = 0
        coll_time = 1e+10
        coll_time_particle = coll_time
        domain_lab = 0
        coll_with_wall = False
        coll_with_particle = False
                        
        #Checking collision time with particles in main domain
        for n in range(0, N):
            for i in range(n+1, N):
                rab = [x[n] - x[i], y[n] - y[i]]
                vab = [v_x[n] - v_x[i], v_y[n] - v_y[i]]
                Disc = np.dot(rab,vab)*np.dot(rab,vab)-np.dot(vab, vab)*(np.dot(rab, rab)-(Ri[i]+Ri[n])*(Ri[i]+Ri[n]))

                if Disc > 0:
                    coll_time_particle = (-np.dot(rab, vab) - sqrt(Disc)) / np.dot(vab, vab)
                else:
                    coll_time_particle = 3e+8

                if coll_time_particle < coll_time and coll_time_particle >= 0:
                    coll_time = coll_time_particle
                    coll_partner_1 = n
                    coll_partner_2 = i
                    coll_with_wall = False
                    coll_with_particle = True
                        
            coll_time_r_wall_x = (E_wall - Ri[n] - x[n]) / (v_x[n] + 1e-20)
            if coll_time_r_wall_x >= 0 and coll_time_r_wall_x <= coll_time:
                coll_time = coll_time_r_wall_x
                coll_partner_1 = n
                coll_partner_2 = r_wall_x
                coll_with_wall = True
                coll_with_particle = False
                
            coll_time_l_wall_x = -(x[n] - Ri[n] - W_wall) / (v_x[n] + 1e-20)
            if coll_time_l_wall_x >= 0 and coll_time_l_wall_x <= coll_time:
                coll_time = coll_time_l_wall_x
                coll_partner_1 = n
                coll_partner_2 = l_wall_x
                coll_with_wall = True
                coll_with_particle = False

            coll_time_u_wall_y = (N_wall - Ri[n] - y[n]) / (v_y[n] + 1e-20)
            if coll_time_u_wall_y >= 0 and coll_time_u_wall_y <= coll_time:
                coll_time = coll_time_u_wall_y
                coll_partner_1 = n
                coll_partner_2 = u_wall_y
                coll_with_wall = True
                coll_with_particle = False
                
            coll_time_l_wall_y = -(y[n] - Ri[n] - S_wall) / (v_y[n] + 1e-20)
            if coll_time_l_wall_y >= 0 and coll_time_l_wall_y <= coll_time:
                coll_time = coll_time_l_wall_y
                coll_partner_1 = n
                coll_partner_2 = l_wall_y
                coll_with_wall = True
                coll_with_particle = False

        #Update positions and velocities 
        if coll_time < dt:
            #Update positions to the point of collision minus some small, negligible,  numerical value to
            #prevent particles from getting stuck in a wall
            x = x + v_x*coll_time*(1-err)
            y = y + v_y*coll_time*(1-err)
          
            #track particle zero
            x_track = x_track + v_x[0]*coll_time*(1-err)
            y_track = y_track + v_y[0]*coll_time*(1-err)

            time_t = time_t + coll_time

            if coll_with_particle == True:
                #Update velocities of colliding particles
                ra = np.array([x[coll_partner_1], y[coll_partner_1]])
                va = np.array([v_x[coll_partner_1], v_y[coll_partner_1]])
                rb = np.array([x[coll_partner_2], y[coll_partner_2]])
                vb = np.array([v_x[coll_partner_2], v_y[coll_partner_2]])
                n = (ra-rb)/sqrt(np.dot((ra-rb), (ra-rb)))
                       
                del_v = n*np.dot((va - vb), n)
                      
                #Update velocities                      
                v_x[coll_partner_1] = v_x[coll_partner_1] - del_v[0]*2*mi[coll_partner_2]/(mi[coll_partner_1]+mi[coll_partner_2])
                v_y[coll_partner_1] = v_y[coll_partner_1] - del_v[1]*2*mi[coll_partner_2]/(mi[coll_partner_1]+mi[coll_partner_2])
                v_x[coll_partner_2] = v_x[coll_partner_2] + del_v[0]*2*mi[coll_partner_1]/(mi[coll_partner_1]+mi[coll_partner_2])
                v_y[coll_partner_2] = v_y[coll_partner_2] + del_v[1]*2*mi[coll_partner_1]/(mi[coll_partner_1]+mi[coll_partner_2])

            elif coll_with_particle == False:
                if coll_partner_2 == r_wall_x or coll_partner_2 == l_wall_x:
                    v_x[coll_partner_1] = -v_x[coll_partner_1]
                if coll_partner_2 == u_wall_y or coll_partner_2 == l_wall_y:
                    v_y[coll_partner_1] = -v_y[coll_partner_1]               

        elif coll_time >= dt: #Update positions and velocities using dt
            x = x + v_x*dt
            y = y + v_y*dt

            #track particle zero
            x_track = x_track + v_x[0]*dt
            y_track = y_track + v_y[0]*dt     
        
            time_t = time_t + dt

        #Check if particles are within bounds
        in_boundaries = True
        for n in range(0, N):
            in_boundaries = (x[n] <= E_wall) and (x[n] >= W_wall) and (y[n] <= N_wall) and (y[n] >= S_wall)
            if not in_boundaries:
                print('outside boundaries')
                sys.exit()
        print('in boundaries')
        print(in_boundaries)
        
        #Check overlap
        there_is_overlap = False
        for n in range(0, N):
            for i in range(0, N):
                if i != n:
                    dx = x[i] - x[n]
                    dy = y[i] - y[n]
                    if dx*dx + dy*dy < (Ri[i]+Ri[n])*(Ri[i]+Ri[n]) - 1e-8:
                        there_is_overlap = True
                    if there_is_overlap == True:
                        print("there is overlap")
                        overlap_counter = overlap_counter + 1
                        sys.exit()
        print(there_is_overlap)
        print('overlap counter')
        print(overlap_counter)

    #Exporting data to file
    my_file = 'BMW_data_norm_dist_' + str(sim) + '.txt'
    print(my_file)
    file = open(my_file, "w")
    for index in range(0, size(t_track)):
        data_r2 = str(r2_track[index])
        data_time = str(t_track[index])
        file.write(data_time)
        file.write(' ')
        file.write(data_r2)
        file.write('\n')

    file.close()
    print("--- %s seconds ---" % (time.time() - start_time))
