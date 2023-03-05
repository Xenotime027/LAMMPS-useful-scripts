#This is a specific script. Use with discretion. The main takeaway from this script is the delienated code (inside the ### section), 
#which introduces the logic behind discrete steps to create a stacking fault (to maintain "good" translational symmetry) along the 112 crystallographic plane.

#This dictionary contains: the keys: alloy name, the items (in order): composition of the first element, number of elements, string of elements
gitpar = [3.6215814614347517, 3.5845818198019868, 3.5965846604277805]
alloy_list = {"MnNi": ["0.6", 2, "Mn,Ni"], "FeCo": ["0.5", 2, "Fe,Co"], "Cantor": ["0.2", 5, "Ni,Fe,Co,Cr,Mn"]}

for num, i in enumerate(alloy_list):
    #Specifies the lattice parameter
    a = gitpar[num]

    #The following lines determine how big your structure file is
    dupx = 6 #repetition along x-axis
    dupy = 9 #repetition along y-axis
    dupz = 8 #repetition along z-axis
    
    #This command creates the probes, not yet alloyed (pure nickel) with the specified lattice constant from before (see above)
    !atomsk --create fcc $a Ni orient [11-2] [-110] [111] -duplicate $dupx $dupy $dupz "$i"_bulk.lmp lmp
    !mv "$i"_bulk.lmp lmp."$i"_bulk
    
    k = alloy_list[i][0] #concentration of the first element (the rest is (1-k)/(n-1))
    n = alloy_list[i][1] #number of elements
    e = alloy_list[i][2] #the list of elements (to be written in LAMMPS file)

    #This line alloys the pure Ni system from before
    #ovitos, my ovitos_code, input file, concentration, out file, tolerance == 0, number of elements, elements list
    
    #Use another code to randomly shuffle the atom types
    !~/LAMMPS/ovitos ~/LAMMPS/ovitos_HEA_types_wrap.py lmp."$i"_bulk $k lmp."$i"_bulk 0 $n $e
    
    #This line inserts vacuum above and below the alloys
    #ovitos, my ovitos_code, input file, vacuum ratio (unitless), out file
    h = 1+0.4
    !~/LAMMPS/ovitos ~/LAMMPS/ovitos_insert_vacuum.py lmp."$i"_bulk $h lmp."$i"_bulk
    
    #################################################################################################################################################
    #alpha_liste: This list specifies how much shift (along the x-axis) is supplied to the stacking fault code. 
    #The shift distance is motivated by crystallographic reasons to maintain appropriate translational symmetry.
    #The factor sqrt 6 corresponds to the 112 crystallographic direction
    alpha_liste = np.linspace(0, 1.5*a/np.sqrt(6), 10) #creates 10 steps along the x-axis direction, multiplied by a calculated distance
    
    for num2, j in enumerate(alpha_liste):
        j = np.round(j, 5) #amount of shift, rounded to 5 decimal places to avoid floating point error
        #Executes the "shifting" code to create stacking fault(s). The shift is directed along x direction on two different heights: 0.33 and 0.67
        !~/LAMMPS/ovitos ~/LAMMPS/ovitos_stapelfehler.py lmp."$i"_bulk lmp."$i"_sf_"$num2" $j
    #################################################################################################################################################
