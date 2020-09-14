import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plab
import scipy.stats
from numpy import *
import matplotlib.animation as manimation
import time

start_time = time.time()

number_of_files = 0#0
start_file = 0
my_file = 'BMW_data_norm_dist_' + str(start_file) + '.txt'
avg_time_data, avg_r2_data = np.genfromtxt(my_file, unpack=True)

for index in range(start_file+1, number_of_files+1):
    my_file = 'BMW_data_norm_dist_' + str(index) + '.txt'
    time_data, r2_data = np.genfromtxt(my_file, unpack=True)
    avg_time_data = avg_time_data + time_data
    avg_r2_data = avg_r2_data + r2_data

avg_time_data = avg_time_data / (number_of_files + 1 - start_file)
avg_r2_data = avg_r2_data / (number_of_files + 1 - start_file)

size_data = size(avg_time_data)

avg_time_data_restricted = avg_time_data[10:size_data]
avg_r2_data_restricted = avg_r2_data[10:size_data]

lin_eq_restr = np.polyfit(avg_time_data_restricted, avg_r2_data_restricted, 1)
p_lin_restr = np.poly1d(lin_eq_restr)

plab.plot(avg_time_data_restricted,p_lin_restr(avg_time_data_restricted),"r--")

result = scipy.stats.linregress(avg_time_data_restricted, avg_r2_data_restricted)
print(result.slope)
print(result.intercept)
print(result.rvalue)
cont = plt.plot(avg_time_data, avg_r2_data)
plt.show()
