import pandas as pd
from math import isnan 
import libchebipy
import re
import argparse


# parse arguments and save the path to the data
parser = argparse.ArgumentParser()
parser.add_argument("dataPath", help="Specify the path to the data file that contains the raw metabolite table")
args = parser.parse_args()
dataPath = args.dataPath

# list to save the index of the rows that have no chebi
emtpyChebis = []

# read data
data = pd.read_csv(dataPath, delimiter="\t")

# process each row in the data
for i in data.index:
    chebi = data['Chebi'][i]

    if (isinstance(chebi, str)):
        # take a single chebi id
        chebiID = re.search('CHEBI:([0-9]+).*', chebi).group(1)

        # look for current chebi id in the ontology
        chebiData = libchebipy.ChebiEntity(chebiID)

        # get all the IDs associated with current chebi
        allChebiIDs = chebiData._ChebiEntity__all_ids

        # check if there is more than one ID associated with the current chebi ID
        if(allChebiIDs is not None):
            # get the first one from the list, which is the main one
            mainID = allChebiIDs[0]

            # update chebi ID value to keep only the main one
            data.loc[[i], ['Chebi']] = "CHEBI:" + str(mainID)
        else:
            data.loc[[i], ['Chebi']] = "CHEBI:" + str(chebiID)
    else:
        # check if instead of chebi there is a "NaN" value
        if(isnan(chebi)):
            # add current index to list of empty chebis 
            emtpyChebis.append(i)

# remove rows without chebi
data.drop(labels = emtpyChebis)

# save processed data in file
newFilePath = dataPath[:len(dataPath)-4] + "_Processed.tsv"
data.to_csv(newFilePath, sep="\t", index=False)


