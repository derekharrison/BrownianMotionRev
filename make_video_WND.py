import numpy as np
import matplotlib.pyplot as plt
from numpy import *
import matplotlib.animation as manimation
import time
import sys

start_time = time.time()

plt.rcParams['animation.ffmpeg_path']='C:/Users/d-w-h/Downloads/ffmpeg-20200818-1c7e55d-win64-static/ffmpeg-20200818-1c7e55d-win64-static/bin/ffmpeg.exe'
writer=manimation.FFMpegWriter(bitrate=20000, fps=15)

fig = plt.figure(figsize=(8,8))

track_x = []
track_y = []

def animate(i):
    track_overlap = []
    overlap_x_vec = []
    overlap_y_vec = []
    my_file_x = 'BMW_norm_dist_x_' + str(i) + '.txt'
    my_file_y = 'BMW_norm_dist_y_' + str(i) + '.txt'
    print(i)
    fig.clear()
    data_x = np.genfromtxt(my_file_x, unpack=True)
    data_y = np.genfromtxt(my_file_y, unpack=True)

    x_zero = data_x[0]
    y_zero = data_y[0]
    track_x.append(x_zero)
    track_y.append(y_zero) 
    ax = plt.axes(xlim=(-3, 3), ylim=(-3, 3))
    cont = plt.scatter(data_x, data_y, s=300, c='green')
    cont = plt.scatter(x_zero, y_zero, s=1500, c='blue')
    cont = plt.plot(track_x, track_y)
    return cont

size_t = 99
anim = manimation.FuncAnimation(fig, animate, frames=size_t, repeat=False)

print("Done Animation, start saving")

anim.save('BMW_norm_dist_yt.mp4', writer=writer, dpi=200)
    
print("--- %s seconds ---" % (time.time() - start_time))
