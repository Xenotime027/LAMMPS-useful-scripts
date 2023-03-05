#A script to create stacking fault (https://en.wikipedia.org/wiki/Stacking_fault) by shifting two atomic planes.
#The shift is given by the user in Angstroms. Requires Ovito.
#Arguments:
#1. Path of the structure file
#2. Path for the export
#3. Value of the shift (Angstrom)

import ovito
from ovito.io import import_file, export_file
from ovito.data import *
from ovito.modifiers import ExpressionSelectionModifier, AffineTransformationModifier, WrapPeriodicImagesModifier
from ovito.data import SimulationCell, DataCollection

import numpy as np
import sys

####################################################################
def shrink_wrap(frame: int, data: DataCollection):
    # There's nothing we can do if there are no input particles.
    if not data.particles or data.particles.count == 0: return

    # Compute min/max range of particle coordinates.
    coords_min = np.amin(data.particles.positions, axis=0)
    coords_max = np.amax(data.particles.positions, axis=0)
    # Build the new 3x4 cell matrix:
    #   (x_max-x_min  0            0            x_min)
    #   (0            y_max-y_min  0            y_min)
    #   (0            0            z_max-z_min  z_min)
    matrix = np.empty((3,4))
    matrix[:,:3] = np.diag(coords_max - coords_min)
    matrix[:, 3] = coords_min

    # Assign the cell matrix - or create whole new SimulationCell object in
    # the DataCollection if there isn't one already.
    SimulationCell.create(data, matrix, (False, False, False))
####################################################################

try:
  input_file = sys.argv[1]
  export_name=sys.argv[2]
  val = sys.argv[3]
except:
  print("Error in the inputs. Please check that all of the inputs must be present: structure file path, export file path, and the value for shifting (Angstrom) ")
  sys.exit(1)

pipeline = import_file(input_file)

#Stacking fault is created by shifting a selected row of atoms by the given value val
#First we shift the above two thirds of the atoms.
expression = f"ReducedPosition.Z>=0.33"
pipeline.modifiers.append(ExpressionSelectionModifier(expression = expression))
pipeline.modifiers.append(AffineTransformationModifier(operate_on={"particles"}, only_selected=True,
                            transformation=[[1,0,0,val], [0,1,0,0], [0,0,1,0]]))

#Then, we shift the rest one third of the atoms.
expression = f"ReducedPosition.Z>=0.66"
pipeline.modifiers.append(ExpressionSelectionModifier(expression = expression))
pipeline.modifiers.append(AffineTransformationModifier(operate_on={"particles"}, only_selected=True,
                            transformation=[[1,0,0,val], [0,1,0,0], [0,0,1,0]]))

#Shrink Wrap to remove vacuum (can be commented out if necessary)
pipeline.modifiers.append(shrink_wrap)

#This line should not be commented out to ensure periodicity.
pipeline.modifiers.append(WrapPeriodicImagesModifier())

export_file(pipeline, export_name, "lammps/data", atom_style="atomic")
