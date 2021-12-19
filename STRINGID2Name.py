import gzip,os

aliasfile = "gz/39947.protein.aliases.v11.5.txt.gz"
interaction_file = "txt/converted_ppi.txt"
output_file = "txt/output_file.txt"
preferred_source = "BLAST_UniProt_GN_Name"
backup_source = "Uniprot"
alias_key = {}

with gzip.open(os.path.join(aliasfile),'rt') as file:
    file.readline()
    line = file.readline()

    while line:
        items = line.rstrip().split('\t')
        tmp = alias_key.get(items[0], {})
        source = items[-1]
        tmp[source] = items[1]
        alias_key[items[0]] = tmp
        line = file.readline()

out = open(os.path.join(output_file),'w')

with open(os.path.join(interaction_file),'r') as file:
    line = file.readline()

    while line:
        items = line.rstrip().split(' ')
        # For each STRING_ID, get the protein name from the preferred source in the alias file
        # If that is missing then go for the backup source
        left_prot = alias_key[items[0]].get(preferred_source,alias_key[items[0]][backup_source])
        right_prot = alias_key[items[1]].get(preferred_source,alias_key[items[1]][backup_source])
        out.writelines(' '.join([left_prot,right_prot,items[-1],'\n']))
        line = file.readline()

out.close()