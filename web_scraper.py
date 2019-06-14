import sys
from Bio import Entrez
from parameter import Parameter

def search(parameter):
    #query to pubmed
    #search and retrieve primary IDs and term translations and optionally retains results
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
        parameter.create_output_directory()
        f=open(parameter.scrapper_result ,"w", encoding='utf-8')
        f.write(data)
        f.close()
        print ("Search results stored in " + parameter.scrapper_result)

def scrape(parameter):
    print("----Step 1: Scraping pubmed database")
    searchResults = search(parameter)
    fetchResults = fetch_details(searchResults['IdList'])
    write_xml(fetchResults, parameter)
    print("----End step 1\n")
    return fetchResults
