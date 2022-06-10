"""

Author: Janith Weeraman
Date: 25/11/2021

Get only the interactions that are above a certain combined score cutoff

"""
import os.path
import gzip

# path to ppi file
path = os.getcwd()
ppi = "gz/39947.protein.links.detailed.v11.5.txt.gz" # ppi file
output = "txt/ppi.tsv"
cutoff = 0.40


newfile = open(os.path.join(path,output),'w')

# Read the original interaction file and write into new tsv file based on cutoff
with gzip.open(os.path.join(path,ppi),'rt') as file:
    # Loop through all lines in a file and get only protein names and combined score
    newcutoff = cutoff * 1000
    while True:
        line = file.readline()
        if line:
            line = line.rstrip().split(' ')
            formatted = [line[i] for i in [0,1,-1]]
            try:
                # Only write files that are above the cutoff
                if int(formatted[-1]) >= newcutoff: newfile.writelines('\t'.join(formatted)+'\n')
            except ValueError:
                # Exception for the first line that contains header data
                newfile.writelines('\t'.join(formatted) + '\n')

        else: break # break out of the loop after last line is written

newfile.close()