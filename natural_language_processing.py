import nltk
from bs4 import BeautifulSoup
from article import Article
from gensim import corpora

"""
Checks nltk modules
"""
def nltk_check():
    #check if punkt and stopwords resources are avaliable
	try:
		nltk.data.find('tokenizers/punkt')
		nltk.data.find('corpora/stopwords')
	except LookupError:
		print("Downloading punkt...")
		nltk.download('punkt')
		print("Downloading stopwords...")
		nltk.download('stopwords')

"""
Parse xml results file to a list of Article objects
Input  : path to xml, result PubMed IDs
Output : Article objects list
"""
def parse_xml(xml_path, pmids):
    article_list = []
	# open and parse xml with BeautifoulSoup
    infile = open(xml_path,"r")
    contents = infile.read()
    soup = BeautifulSoup(contents,'xml')
	try:
		# parse xml file into a list of Article objects
		titles = soup.find_all('ArticleTitle')
		abstracts = soup.find_all('Abstract')
		print("* Parsing xml file...")
		for i in range(0, len(titles)):
			article_list.append(Article(pmids[i], titles[i].text, abstracts[i].text))
	except:
		print("An error ocurred while parsing xml file. Try it again")
		sys.exit(1)
    return article_list

def preprocessing(article_list):
    doc_arrays = []
    for article in article_list:
        doc_arrays.append(article.get_abstract_array())

    return doc_arrays

def process_docs(parameter, pmids):
    #Check if necessary resources are avaliable
    nltk_check()

    #Parse xml file
    article_list = parse_xml(parameter.data_extraction_result, pmids)
    doc_arrays = preprocessing(article_list)

    #Creating dictionary and corpus
    dictionary = corpora.Dictionary(doc_arrays)
    dictionary.save(parameter.dictionary)

    corpus = [dictionary.doc2bow(array) for array in doc_arrays]
    corpora.MmCorpus.serialize(parameter.corpus, corpus)

    return article_list
