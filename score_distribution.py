import os.path

# path to ppi file
path = "C:\\Users\\Janith\\Documents\\"
ppi = "39947.protein.links.full.v11.5.txt" # ppi file


with open(os.path.join(path,ppi),'r') as file:
    file = [line.rstrip() for line in file.readlines()]

#titles = file[0].split(' ')[2:]
#data = [line[2:] for line in file[1:]]

for line in file:
    line = line.split(' ')
    if line[9] > 0: print(line)