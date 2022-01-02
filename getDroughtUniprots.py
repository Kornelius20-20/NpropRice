changed = open("txt/processed_uniprot.txt", 'w')

with open("txt/processed_uniprot.csv", 'r') as file:
    # Copy title line to file
    line = file.readline()
    items = line.rstrip().split(',')
    changed.writelines('\t'.join(items) + '\n')

    # start at next line
    # line = file.readline()
    # For each line
    while line:
        items = line.rstrip().split('\t')
        # items = items[2:]
        # items[0] = items[0].split(',')[0]
        if "drought" in line.lower(): changed.writelines('\t'.join(items) + '\n')
        line = file.readline()

changed.close()
