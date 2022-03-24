def rescale(valuedict,scalemax=256,mincutoff=0):
    """
    Method to rescale a dictionary with numeric values

    :param valuedict: the dictionary with keys and numeric values
    :param scalemax: the max value to scale to
    :param mincutoff: the value below which inputs will be rounded to 0
    :return:
    """

    scalemax = scalemax -1
    smol = min(valuedict.values())
    big = max(valuedict.values())
    for key in valuedict.keys():
        newval = int(valuedict[key]/(big-smol) * scalemax)
        valuedict[key] = newval if newval >= mincutoff else 0

"""

Using the GO REST API to get the related GO terms for a particular term

"""

def get_go_json(goterm):
    import requests, sys, json
    requestURL = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{goterm}/descendants?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates"

    r = requests.get(requestURL, headers={"Accept": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    responseBody = r.json()

    with open('response.json', 'w') as file:
        json.dump(responseBody, file)


if __name__ == "__main__":
    goterm = "GO:0009819" # drought response
    filename = "response.json"
    get_go_json(goterm)