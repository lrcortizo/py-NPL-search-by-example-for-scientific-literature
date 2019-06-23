import sys
from Bio import Entrez
from parameter import Parameter

def search(parameter):
    #query to pubmed
    #search and retrieve primary IDs and term translations and optionally retains results
    Entrez.email = 'your.email@example.com'
    try:
        searchHandle = Entrez.esearch(db='pubmed',
                                retmax=parameter.max_results,
                                retmode='xml',
                                term=parameter.search_term)
        searchResults = Entrez.read(searchHandle)
        searchHandle.close()
        return searchResults
    except:
        return None

def fetch_details(id_list):
    #retrieve records in the requested format from a list of IDs
    ids = ','.join(id_list)
    Entrez.email = 'your.email@example.com'
    try:
        fetchHandle = Entrez.efetch(db='pubmed',
                               retmode='xml',
                               id=ids)
        fetchResults = fetchHandle.read()
        fetchHandle.close()
        return fetchResults
    except:
        return None

def write_xml(data, parameter):
    #check and write results
    if data==None:
        print (80*"*"+"\n")
        print ("This search returned no hits")
        sys.exit(1)

    else:
        #write result into a xml file
        f=open(parameter.data_extraction_result ,"w", encoding='utf-8')
        f.write(data)
        f.close()
        print ("Search results stored in " + parameter.data_extraction_result)

def extract(parameter):
    searchResults = search(parameter)
    fetchResults = fetch_details(searchResults['IdList'])
    write_xml(fetchResults, parameter)

    return fetchResults
