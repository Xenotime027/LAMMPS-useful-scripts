#This code is very specific to my problem
#It depends on the parameters of the LAMMPS simulations you used, such as the timestep value
#Please use this if you know what you are doing
#Calculates the mean square difference from a LAMMPS trajectory file from each timestep 
#by comparing the atomic positions across two different frames

import ovito
from ovito.io import import_file, export_file
from ovito.data import *
import numpy as np
import pandas as pd
import sys

import ase
from ase.data import atomic_names, atomic_masses, vdw_radii
from ovito.modifiers import CalculateDisplacementsModifier as CDM, ExpressionSelectionModifier as ESM, DeleteSelectedModifier as DSM

print("Hello, this is OVITO %i.%i.%i" % ovito.version)
input_file = str(sys.argv[1])
export_name = str(sys.argv[2])

try:
    frame_interval = int(sys.argv[3]) #sets your interval, how often you calculate the msd
except:
    frame_interval = 1 #default value for interval, skips no frame

n_elem = np.arange(5) #number of elements, please change if your file has different values
emesde_all = []

for num in n_elem: #iterates through all atom types
    pipeline = import_file(input_file, sort_particles=True)
    my_expression = f"ParticleType != {num+1}"
    pipeline.modifiers.append(ESM(expression=my_expression))
    pipeline.modifiers.append(DSM())
    data = pipeline.compute()
    n_atoms = data.particles.count

    modifier = CDM(use_frame_offset=True, frame_offset=-1)
    pipeline.modifiers.append(modifier)
    length = len(range(0,pipeline.source.num_frames,frame_interval))
    length = length - 100 #this depends on your simulation parameters!

    emesde = np.zeros((length, n_atoms))
    for frame in range(101, pipeline.source.num_frames, frame_interval):
        data = pipeline.compute(frame)
        emesde[frame-100, :] = emesde[frame-100-1, :] + (data.particles["Displacement Magnitude"]**2)

    emesde_sum = np.sum(emesde, axis=1)/len(data.particles["Displacement Magnitude"])
    emesde_all.append(emesde_sum)

np.savetxt(export_name, emesde_all)
