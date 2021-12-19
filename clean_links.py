"""

Author: Janith Weeraman
Date: 25/11/2021

"""
import os.path
import gzip

# path to ppi file
path = os.getcwd()
ppi = "gz/39947.protein.links.detailed.v11.5.txt.gz" # ppi file

cutoff = 0.8

newfile = open(os.path.join(path,"txt/new_ppi.txt"),'w')

with gzip.open(os.path.join(path,ppi),'rt') as file:
    # Loop through all lines in a file and get only protein names, experimental score and combined score
    while True:
        line = file.readline()
        if line:
            line = line.rstrip().split(' ')
            formatted = [line[i] for i in [0,1,-1]]
            try:
                if int(formatted[-1])/1000 >= cutoff: newfile.writelines(' '.join(formatted)+'\n')
            except ValueError:
                newfile.writelines(' '.join(formatted) + '\n')

        else: break

newfile.close()