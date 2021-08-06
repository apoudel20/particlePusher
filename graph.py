# import ast # parse list from file

# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# from mpl_toolkits import mplot3d

# fig = plt.figure()
# ax = plt.axes(projection='3d')
# # ax.set_xlim([0,100])
# # ax.set_ylim([0,100])
# # ax.set_zlim([0,100])

# for i in range(6):

#     position_data = pd.read_csv("frames/output_"+str(i)+".csv")
    
#     x = position_data.iloc[:,0][::1000]    
#     y = position_data.iloc[:,1][::1000]
#     z = position_data.iloc[:,2][::1000]
#     ax.scatter3D(x, y, z, c = z)

# plt.show()

import ast # parse list from file

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# from numba import vectorize

frames_rendered = 16
x,y,z = [],[],[]

for i in range(frames_rendered):
    position_data = pd.read_csv("frames/output_"+str(i)+".csv")

    x.append(position_data.iloc[:,0][::1000])   
    y.append(position_data.iloc[:,1][::1000])
    z.append(position_data.iloc[:,2][::1000])


def animate_vector(i):
    title.set_text("Frame: " + str(i))
    position_plots._offsets3d = (x[i], y[i], z[i])


fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

ax.set_xlim([-15,15])
ax.set_ylim([-5,15])
ax.set_zlim([-15,15])
title = ax.set_title('3D Test')
def get_cube():   
    phi = np.arange(1,10,2)*np.pi/4
    Phi, Theta = np.meshgrid(phi, phi)

    x = np.cos(Phi)*np.sin(Theta)
    y = np.sin(Phi)*np.sin(Theta)
    z = np.cos(Theta)/np.sqrt(2)
    return x,y,z


a = 25
b = 10
c = 25
cube_x,cube_y,cube_z = get_cube()


position_plots = ax.scatter3D(x[0],y[0],z[0])
        
myAnimation = FuncAnimation(fig, animate_vector, frames=np.arange(1,frames_rendered), blit = False, interval=1)
ax.plot_surface(cube_x*a+1, cube_y*b+b/2, cube_z*c+1,alpha=0.1)

plt.show()

