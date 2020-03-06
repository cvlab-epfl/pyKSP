
import sys
sys.path.insert(0, 'pyKSP')

import numpy as np
import matplotlib.pyplot as plt


# Define a couple of usefull functions
def get_fake():
    # To create fake tracks
    out = np.zeros((37,36)) + 0.00001
    out[10,10] = 0.99
    out[20,20] = 0.99
    
    return out

def convert_v(v,dims,depth):
    # To convert the output of KSP into tracks format
    MAP = np.zeros(dims)
    n_tracks = len(v) // depth
    out = np.zeros((depth,n_tracks,2))
    for t in range(depth-1):
        for i in range(n_tracks):
            out[t,i,0] = v[t*n_tracks + i]/dims[1]
            out[t,i,1] = v[t*n_tracks + i]%dims[1]
            
    return np.int32(out)


# Import the wrapper
from ksp import pyKShorthestPathGraph

'''
    We will need to create a pyKShorthestPathGraph object. The constructor takes the following inputs, which we are going to build:
    - A flattened vector containing the costs of going through each grid location at each time step.
    - The W and H dimensions of your grid here 36 and 37.
    - The number of frames to be processed, called depth, here 10.
    - The admissible radius for two consecutive detections of the same track.
    - The access points of your grid, in flattened format as well. These are the points where tracks can enter and exit the grid.
'''

depth = 30 #Number of frames

# Prepare flat vector for ksp containing costs at each grid locations at each time step
q_vector = np.zeros(37*36*depth)
for t in range(0,depth):
    Q_loc = get_fake()
    flat_q = np.clip(np.ndarray.flatten(Q_loc),1e-6,0.999999)
    q_vector[37*36*t:37*36*(t+1)] = -np.log(flat_q/(1-flat_q)) # Costs in -log() format
    
access_points = np.asarray([0]) # Define the access points on your grid
G = pyKShorthestPathGraph(q_vector,36,37,depth,4,access_points) #Be carefull ordering of dimensions inverted
v = G.getPath(0,depth-1) # From_frame - To_frame (inclusive). Be carreful if you set To_frame>depth - 1, you get a memory leak
print(v)
del G


# Convert output to a "tracks" format
out = convert_v(v,(37,36),depth)
print(out)

