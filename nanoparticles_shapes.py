import ase
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ase.cluster import Icosahedron as ico, Decahedron as deca, Octahedron as octa

def create_ico(N, symbol="Ni", noshells=8, latticeconstant=3.524):
    #Creates an icosahedral nanoparticle with a given number of shells and lattice constant
    structure = ico(symbol, noshells, latticeconstant)
    ase.io.write(f"lmp.1_{N}", structure, format="lammps-data")
    
def create_octa(N, symbol="Ni",n=2, latticeconstant=3.524):
    #Creates a regularly truncated octahedral nanoparticle with a given number of shells and lattice constant
    length = n*3+1
    cutoff = (length-1)//3
    structure = octa(symbol, length, cutoff, latticeconstant)
    ase.io.write(f"lmp.2_{N}", structure, format="lammps-data")

def create_octa_nt(N, symbol="Ni",n=2, latticeconstant=3.524):
    #Creates a non-truncated octahedral nanoparticle with a given number of shells and lattice constant
    length = n*3+1
    cutoff = 0
    structure = octa(symbol, length, cutoff, latticeconstant)
    ase.io.write(f"lmp.3_{N}", structure, format="lammps-data")
    
def create_deca(N, symbol="Ni", p=3, q=3, r=3, latticeconstant=3.524):
    #Creates a Marks decahedral nanoparticle with set symmetry parameters (3 values must be set) and lattice constant
    structure = deca(symbol, p, q, r, latticeconstant)
    ase.io.write(f"lmp.4_{N}", structure, format="lammps-data")
    
def create_deca_nt(N, symbol="Ni", p=3, q=3, r=3, latticeconstant=3.524):
    #Creates a non-truncated decahedral nanoparticle with set symmetry parameters (3 values must be set) and lattice constant
    structure = deca(symbol, p, q, r, latticeconstant)
    ase.io.write(f"lmp.5_{N}", structure, format="lammps-data")
