#This script allows changing and adding element types in a monoatomic LAMMPS structure file. Requires Atomic Simulation Environment (ASE) library.
import numpy as np
from ase.data import *
import sys
import os

####################################################################################################
def overwrite_type_lammps(filename, elements_list):
  #This function changes the type of a given LAMMPS structure file by adding the types determined by the user.

  #Input-Parameters:
  #filename: file path of the structure file to be edited.
  #elements_list: placeholder for the index, element symbol, and atomic masses.

  #Extract the content
  number, element, mass = zip(*elements_list)
  print(f"Number of elements: {len(number)}")
  
  #Begin editing the files
  with open(filename, "r+") as f:
    c = f.read().find("atom types")
    f.seek(c-2)
    #Edits the number of atom types
    f.write(f"{len(number)} atom types\n")

    f.seek(0)
    c = f.read().find("Atoms  # atomic")
    f.seek(c)
    old = f.read()

    f.seek(0)
    c = f.read().find("Masses")
    f.seek(c)
    #Edits the details: masses and element symbols
    f.write("Masses\n\n")
    for i in np.arange(len(number)):
      f.write(f"{i+1} {mass[i]}  # {element[i]}\n")
      f.write("\n"+old)
  f.close()
####################################################################################################    
   
#Implementation
input_file = sys.argv[1]
temp_file = input_file+"_temp"
os.system(f"cp {input_file} {temp_file}")
elements_list = []

#Define here the element list to be added, can be customised to accept inputs from the terminal by sys.argv()
elements = ["Ni", "Fe", "Co", "Cr", "Mn"]

#Creates the list of elements with index, symbol, and masses.
for i in np.arange(len(elements)):
  elements_list.append([i+1, elements[i], atomic_masses[atomic_numbers[elements[i]]]])

print(elements_list)
print("\nBegin editing the LAMMPS structure file.\n")
overwrite_type_lammps(temp_file, elements_list)
