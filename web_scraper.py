import os
from Bio import Entrez

def search(query):
    Entrez.email = 'your.email@example.com'
    try:
        searchHandle = Entrez.esearch(db='pubmed',
                                sort='relevance',
                                retmax='10',
                                retmode='xml',
                                term=query)
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

def write_xml(data):
    if data==None:
        print (80*"*"+"\n")
        print ("This search returned no hits")

    else:
        if not os.path.exists("tmp/"):
            os.mkdir("tmp/")
        f=open("tmp/web_scrapper_results.xml" ,"w")
        f.write(data)
        f.close()

if __name__ == '__main__':
    searchResults = search('Clenbuterol')
    id_list = searchResults['IdList']
    fetchResults = fetch_details(id_list)
    write_xml(fetchResults)
