import os
import pandas as pd
from suds.client import Client

outdir = "outputs/results"
url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl"

client = Client(url)

client.wsdl.services[0].setlocation(
    'https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')

# authenticate user email
client.service.authenticate('2017s16486@stu.cmb.ac.lk')

idType = 'UNIPROT_ID'
listType = 0

# getGeneClusterReport
overlap = 4
initialSeed = 4
finalSeed = 4
linkage = 0.35
kappa = 50

for _, _, files in os.walk(outdir):
    # Get only text files in the output results directory
    for file in files:
        if 'tsv' in file:
            # with open(os.path.join(outdir,file),'r') as clustfile:
            #     data = [line.rstrip() for line in clustfile.readlines()]
            #     print(data)

            fr = pd.read_csv(os.path.join(outdir, file), delimiter='\t', index_col=False)

stuff = []
for col in fr.columns:
    listName = col
    out = fr[col].dropna().tolist()
    for i in out: stuff.append(i)
# inputIds = ",".join(out)

inputIds = ",".join(set(stuff))

print(client.service.addList(inputIds, idType, listName, listType))

geneClusterReport = client.service.getGeneClusterReport(overlap, initialSeed, finalSeed, linkage, kappa)

totalClusters = len(geneClusterReport)
print('Total groups:', totalClusters)
resF = 'geneClusterReport.txt'
with open(resF, 'w') as fOut:
    for simpleGeneClusterRecord in geneClusterReport:
        # EnrichmentScore = simpleGeneClusterRecord.score
        fOut.write(simpleGeneClusterRecord.name + '\tEnrichmentScore: ' + str(simpleGeneClusterRecord.score) + '\n')
        fOut.write(idType + '\tGene Name\n')
        for listRecord in simpleGeneClusterRecord.listRecords:
            gene_id = ','.join(listRecord.values)

            fOut.write(gene_id + '\t' + listRecord.name + '\n')
print('write file:', resF, 'finished!')
