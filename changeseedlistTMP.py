import os,gzip
aliases = "gz/39947.protein.aliases.v11.5.txt.gz"
seedfile = "txt/newseeds.txt"

alias_key = {}
with gzip.open(os.path.join(aliases), 'rt') as file:
    line = file.readline()

    # Create a dict with Alias: Corresponding ID
    while line:
        items = line.rstrip().split()[:2]
        alias_key[items[1]] = items[0]
        line = file.readline()

# Convert seed list to STRING ID list
seedlist = []
with open(os.path.join(seedfile),'r') as file:
    seed = file.readline().rstrip()
    while seed:
        if seed in alias_key.keys(): seedlist.append(alias_key[seed])
        seed = file.readline().rstrip()

print(seedlist)

with open(os.path.join("txt/seed_id.txt"),'w') as file:
    file.writelines('\n'.join(seedlist))
