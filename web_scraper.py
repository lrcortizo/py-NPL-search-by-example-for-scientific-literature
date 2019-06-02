import sys
from Bio import Entrez
from parameter import Parameter

def search(parameter):
    Entrez.email = 'your.email@example.com'
    try:
        searchHandle = Entrez.esearch(db='pubmed',
                                sort='relevance',
                                retmax=parameter.max_results,
                                retmode='xml',
                                term=parameter.search_term)
        searchResults = Entrez.read(searchHandle)
        searchHandle.close()
        return searchResults
    except:
        return None

def fetch_details(id_list):
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
    if data==None:
        print (80*"*"+"\n")
        print ("This search returned no hits")
        sys.exit(1)

    else:
        parameter.create_output_directory()
        f=open(parameter.scrapper_result ,"w")
        f.write(data)
        f.close()
        print ("Search results stored in tmp/web_scrapper_results.xml")

def web_scrapper(parameter):
    searchResults = search(parameter)
    fetchResults = fetch_details(searchResults['IdList'])
    write_xml(fetchResults, parameter)
    return fetchResults
