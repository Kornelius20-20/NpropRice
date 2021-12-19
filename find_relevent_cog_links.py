import gzip

relevent_cogs = "available.txt"
cog_links = "COG.links.v11.5.txt.gz"

with open(relevent_cogs,'r') as file:
    cogs = [line.rstrip().split('\t')[1] for line in file.readlines()]

cogs = set(cogs)
drought_cogs = open('drought_cog_links.txt','w')

with gzip.open(cog_links,'rt') as file:
    line = True
    while line:
        line = file.readline()
        items = line.rstrip().split(' ')
        try:
            if items[0] in cogs or items[1] in cogs:
                drought_cogs.writelines(line)
        except IndexError:
            print(line)
            drought_cogs.close()

drought_cogs.close()
