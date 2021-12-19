
changed = open("processed_uniprot.txt",'w')

with open("uniprot_proteins.txt",'r') as file:
    line = file.readline()
    items = line.rstrip().split('\t')[2:]
    changed.writelines('\t'.join(items) + '\n')
    while line:
        items = line.rstrip().split('\t')
        items = items[2:]
        items[0] = items[0].split(' ')[0]
        if "drought" in line.lower(): changed.writelines('\t'.join(items) + '\n')
        line = file.readline()

changed.close()