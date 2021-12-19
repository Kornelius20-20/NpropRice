import gzip

present_cogs = open('available.txt','w')

with gzip.open("COG.mappings.v11.5.txt.gz",'rt') as file:
    line = True
    i = 0
    while line:
        line = file.readline()
        data = line.rstrip().split("\t")

        if not data[-1] == 'annotation not available' and 'drought' in data[-1].lower():
            present_cogs.writelines('\t'.join((data[0],data[3],data[4])) + '\n')

present_cogs.close()