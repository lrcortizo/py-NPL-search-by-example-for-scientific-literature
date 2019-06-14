import nltk
from bs4 import BeautifulSoup
from article import Article
from gensim import corpora

def nltk_check():
    #check if punkt and stopwords resources are avaliable
	try:
		nltk.data.find('tokenizers/punkt')
		nltk.data.find('corpora/stopwords')
	except LookupError:
		nltk.download('punkt')
		nltk.download('stopwords')

def parse_xml(xml_path):
    #parse xml file into a list of Article objects
    article_list = []

    infile = open(xml_path,"r")
    contents = infile.read()
    soup = BeautifulSoup(contents,'xml')

    pmids = soup.find_all('PMID')
    titles = soup.find_all('ArticleTitle')
    abstracts = soup.find_all('Abstract')

    for i in range(0, len(pmids)):
        article_list.append(Article(pmids[i].text, titles[i].text, abstracts[i].text))
    return article_list

def preprocessing(article_list):
    doc_arrays = []
    for article in article_list:
        doc_arrays.append(article.get_abstract_array())

    return doc_arrays

def process_docs(parameter):
    #Check if necessary resources are avaliable
    nltk_check()

    #Parse xml file
    article_list = parse_xml(parameter.scrapper_result)
    doc_arrays = preprocessing(article_list)

    #Creating dictionary and corpus
    dictionary = corpora.Dictionary(doc_arrays)
    dictionary.save(parameter.dictionary)

    corpus = [dictionary.doc2bow(array) for array in doc_arrays]
    corpora.MmCorpus.serialize(parameter.corpus, corpus)

    return article_list
