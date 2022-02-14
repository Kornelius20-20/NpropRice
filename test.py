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