import pandas as pd
import os,gzip

filepath = r"gz\uniprot-filtered-organism Oryza+sativa+subsp.+japonica+(Rice)+[39947%--.tab.gz"


frame = pd.read_csv(gzip.open(os.path.join(os.getcwd(),filepath)),delimiter='\t')

frame["Annotation"] = [nums[0] for nums in frame["Annotation"]]


# for line in frame["Function [CC]"]:
    # if "stress" in str(line).lower():
    #     print(line)

# Get only rows where the function is known to be related to stress
testframe = frame[frame['Function [CC]'].str.contains('stress',na=False)]
testframe = testframe[testframe['Status'] == 'reviewed']
testframe.to_csv("txt/processed_uniprot.csv")