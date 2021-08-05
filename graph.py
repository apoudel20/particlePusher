import ast # parse list from file

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

fig = plt.figure()
ax = plt.axes(projection='3d')
# ax.set_xlim([0,100])
# ax.set_ylim([0,100])
# ax.set_zlim([0,100])

for i in range(3):

    position_data = pd.read_csv("frames/output_"+str(i)+".csv")
    
    x = position_data.iloc[:,0][::1000]    
    y = position_data.iloc[:,1][::1000]
    z = position_data.iloc[:,2][::1000]
    ax.scatter3D(x, y, z, c = z)

plt.show()
