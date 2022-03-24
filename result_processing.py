import os
import pandas as pd

outputdir = "outputs/results"

def concat_outputs(outputdir):
    outputprots = []
    for dir,folders,files in os.walk(outputdir):
        for file in files:
            fr = pd.read_csv(os.path.join(outputdir,file))
            fr = fr.drop(columns='Unnamed: 0')
            outputprots.append(fr)

    outputprots = pd.concat(outputprots).drop_duplicates('Entry')

    return outputprots


def get_best_csv_name(candset,outputdir):

    filescores = {}
    for dir,folders,files in os.walk(outputdir):
        for file in files:
            with open(os.path.join(outputdir,file),'r') as res:
                # Get only proteins from each csv
                lines = [line.rstrip().split(',')[1] for line in res.readlines()[1:]]

                filescores[file] = len(candset.intersection(set(lines)))

    # return key corresponding to the maximum score
    keys,vals = [],[]
    for key,value in filescores.items():
        keys.append(key)
        vals.append(value)

    bestscore = max(vals)
    outfiles = []
    while True:
        try:
            targetindex = vals.index(bestscore)
            outfiles.append(keys[targetindex])
            vals.pop(targetindex)
            keys.pop(targetindex)
        except ValueError:
            break

    return outfiles


if __name__ == "__main__":
    outputdir = "outputs/results"
    outputprots = concat_outputs(outputdir)
    candset = set(outputprots['Entry'].tolist())
    print(get_best_csv_name(candset,outputdir))