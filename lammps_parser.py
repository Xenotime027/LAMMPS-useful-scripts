#Read / parses the LAMMPS log file.
class lammps_data:
    def __init__(self, file):
        import pandas as pd
        self.file = str(file)
        self.start, self.end, self.atom = [],[],0  #initialises: start, end = row numbers to start and end parsing, atom = initial number of atoms
        self.value = []                            #value = the values to be inserted to the dataframe in self.data
        
        with open(self.file) as file:
            self.text_to_parse = file.readlines()               #sves the read data into the variable text_to_parse
        
        for i in range(len(self.file_to_parse)):                #finds how many different log data files are appended in the given file
            if self.file_to_parse[i].find("Per MPI rank") != -1:  
                self.start.append(i+1)                
            elif self.file_to_parse[i].find("Loop time") != -1:
                self.end.append(i)
            elif self.file_to_parse[i].find("reading atoms ...") !=-1:
                self.atom = self.file_to_parse[i+1].split()[0]

        self.columns = self.file_to_parse[self.start[0]].split()  #parses the column names, must be separated from the numerical contents
        self.start[0] = self.start[0]-1
        
        for i in range(len(self.start)):                      #appends different log files to a placeholder array
            self.cut_file = self.file_to_parse[self.start[i]+2:self.end[i]]
            for line in range(len(self.cut_file)):
                self.value.append(self.cut_file[line].split())  

        self.data = pd.DataFrame(self.value, dtype="float")     #converts into a Pandas Dataframe 
        self.data.columns = self.columns
