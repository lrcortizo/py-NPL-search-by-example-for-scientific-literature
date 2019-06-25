import sys
from Bio import Entrez
from parameter import Parameter

"""
Searchs and get PMIDs from PubMed database
Input  : parameter object
Output : PubMed IDs list
"""
def search(parameter):
    Entrez.email = 'email@example.com'
    print ("* Querying PubMed ...")
    # searchs and retrieves primary IDs, term translations and optionally retains results
    try:
        searchHandle = Entrez.esearch(db='pubmed',
                                retmax=parameter.max_results,
                                retmode='xml',
                                term=parameter.search_term)
        searchResults = Entrez.read(searchHandle)
        searchHandle.close()
    except:
        print ('\033[91m'+20*"*"+" Error occurred while searching PubMed"+'\033[0m')
        sys.exit(1)

    if len(searchResults['IdList']) == 0:
        print ('\033[91m'+20*"*"+" This search returned no hits\n"+'\033[0m')
        sys.exit(1)

    print ("** "+str(len(searchResults['IdList'])) + " articles founded for '"+parameter.search_term+"'")
    return searchResults['IdList']

"""
Retrives the articles from the PubMed IDs
Input  : PubMed IDs
Output : result file in xml format
"""
def fetch_details(pubmedIDs):
    ids = ','.join(pubmedIDs)
    Entrez.email = 'email@example.com'
    try:
        print ("*** Retrieving data ...")
        # retrieves records in the requested format from a list of IDs
        fetchHandle = Entrez.efetch(db='pubmed',
                               retmode='xml',
                               id=ids)
        fetchResults = fetchHandle.read()
        fetchHandle.close()
    except:
        print ('\033[91m'+20*"*"+" Error occurred while retrieving data"+'\033[0m')
        sys.exit(1)
    return fetchResults


"""
Writes the result in a xml file
Input  : Results of querying PubMed
         parameter object
"""
def write_xml(data, parameter):
    try:
        print ("**** Saving data ...")
        #write result into a xml file
        f=open(parameter.data_extraction_result ,"w", encoding='utf-8')
        f.write(data)
        f.close()
        print ("***** Search results stored in " + parameter.data_extraction_result)
    except:
        print('\033[91m'+20*"*"+" The xml file could not be saved"+'\033[0m')
        sys.exit(1)

"""
Main method of extrac_data
Input  : parameter object
Output : Resulting PubMed IDs
"""
def extract(parameter):
    pubmedIDs = search(parameter)
    fetchResults = fetch_details(pubmedIDs)
    write_xml(fetchResults, parameter)

    return pubmedIDs, fetchResults
