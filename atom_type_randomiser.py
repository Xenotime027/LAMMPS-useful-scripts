#This script creates randomised alloy (random solid solution) for a given LAMMPS data file. Requires Ovito and ASE libraries. 
#The user sets the elements and their concentration values, which will be normalised.
#This script takes multiple inputs
#1. Input file location
#2. Export file location
#3. Tolerance rate = 0: gives a buffer on the simulation box (>0 for simulating particles suspended in vacuum)
#4. Optional: name of elements, default Ni-Fe-Co-Cr-Mn
#5. Optional: concentration of each element, default: equiatomic (1/n_elem)

class MemberError(Exception):
    pass
  
import ovito
from ovito.io import import_file, export_file
from ovito.data import *
import numpy as np
import sys

import ase
from ase.data import atomic_names, atomic_masses, vdw_radii, atomic_numbers

input_file = str(sys.argv[1])
export_name = str(sys.argv[2])
try:
    tol_rate = float(sys.argv[3])
except: 
    tol_rate = 0
    print(f"Default tolerance rate is used {tol_rate}.\n")

try:
    names = np.array(sys.argv[4].split(" "))
    n_elem = len(names)
    concentration = np.array(sys.argv[5].split(" "), dtype="float")
    n_conc = len(concentration)
	
    if n_elem != n_conc:
        raise MemberError
        
    concentration = concentration/np.sum(concentration)
    print(f"Normalised concentration: {concentration}")

except MemberError:
    print("Number of elements and concentrations are not same.\n")
    sys.exit(1)

    
except ValueError:
    print("At least one of the concentration values is not numeric.\n")
    sys.exit(1)

except:
    names = ["Ni", "Fe", "Co", "Cr", "Mn"]
    n_elem = len(names)
    concentration = (1/n_elem)*np.ones(n_elem)
    print(f"Default settings:\nNumber of elements: {n_elem}\nElements:{names} \nConcentration: {concentration}")

conc_list = np.array(concentration)

print(conc_list)

pipeline = import_file(input_file)

pipeline.modifiers.append(modify_type)
pipeline.compute()
pipeline.modifiers.append(modify)

if tol_rate!=0:
    pipeline.modifiers.append(modify_wrap)
    
export_name = f"./{export_name}"
export_file(pipeline, export_name, "lammps/data", atom_style="atomic")

########################################
def _create_initial_conc_vec(length, conc):
    conc_vec = []
    for i, c in enumerate(conc):
        conc_vec += ([i+1]*int(length*c))
    return conc_vec

def _pad_conc_vec(conc_vec, length, conc):
    x_p = np.cumsum(conc)
    while len(conc_vec) < length:
        conc_vec.append(np.sum(x_p < np.random.random())+1)
    return conc_vec

def _calculate_concentration(data):
    pt, count = np.unique(data.particles['Particle Type'], return_counts=True)
    count = count / np.sum(count)
    print('New concentrations:')
    for p, c in zip(pt, count):
        print(f'Type {p}: {c}')
    return pt, count

def modify(frame: int, data: DataCollection, conc=(0.5, 0.5), seed=123456,
           only_selected=False):
    assert (not only_selected) | ('Selection' in data.particles.keys()), \
        'No selection defined!'
    assert np.isclose(np.sum(conc), 1), f'sum conc = {np.sum(conc)} != 1'

    # new np syntax (currently not supported by ovito)
    # rng = np.random.default_rng(seed)
    seed = np.random.randint(0, 10000)
    np.random.seed = seed
    conc = conc_list
    print(conc)

    if only_selected:
        select = np.array(data.particles['Selection'])
        select = select.astype(int)

    length = data.particles.count if not only_selected else np.sum(select)
    conc_vec = _create_initial_conc_vec(length, conc)
    conc_vec = _pad_conc_vec(conc_vec, length, conc)

    conc_vec = np.array(conc_vec)
    np.random.shuffle(conc_vec)

    if only_selected:
        data.particles_['Particle Type_'][np.where(select == 1)[0]] = conc_vec
    else:
        data.particles_['Particle Type_'][...] = conc_vec

    _ = _calculate_concentration(data)

from ase.data import atomic_names, atomic_masses, vdw_radii

def modify_type(frame: int, data: DataCollection):
    ptypes =  data.particles_.particle_types_
    masses, radii = [], []
    
    for i in np.arange(len(names)):
        print(i)
        masses.append(atomic_masses[atomic_numbers[names[i]]])
        radii.append(vdw_radii[atomic_numbers[names[i]]])
    
    num_types = len(ptypes.types)
    if num_types < len(names):
        for i in np.arange(num_types, len(names)):
            ptypes.types.append(ParticleType(id = i+1))
            
    for i in np.arange(len(names)):
        type = ptypes.types[i]
        type = ptypes.make_mutable(type)
        type.name = names[i]
        type.mass = masses[i]
        type.vdw_radius = vdw_radii[i]
    
def modify_wrap(frame: int, data: DataCollection):
    # There's nothing we can do if there are no input particles.
    if not data.particles or data.particles.count == 0: return

    # Compute min/max range of particle coordinates.
    coords_min = np.amin(data.particles.positions, axis=0)
    coords_max = np.amax(data.particles.positions, axis=0)

    matrix = np.empty((3,4))
    toleranz = tol_rate*(np.max(coords_max)-np.min(coords_min))
    print(toleranz)
    matrix[:,:3] = np.diag(coords_max - coords_min + toleranz)
    matrix[:, 3] = coords_min - toleranz/2

    # Assign the cell matrix - or create whole new SimulationCell object in
    # the DataCollection if there isn't one already.
    SimulationCell.create(data, matrix, (False, False, False))

#######################################
