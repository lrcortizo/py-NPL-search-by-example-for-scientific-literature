import logging
import nltk
from bs4 import BeautifulSoup
from article import Article
from gensim import corpora
from smart_open import smart_open

def nltk_check():
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
    texts = []
    for article in article_list:
        texts.append(article.get_abstract_array())
    #print(texts)
    return texts

def process_text(parameter):
    print("----Step 2: Text processing")
    #logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
	
    nltk_check()
	
    article_list = parse_xml(parameter.scrapper_result)
    texts = preprocessing(article_list)
    #print(texts)

    dictionary = corpora.Dictionary(texts)
    dictionary.save(parameter.dictionary)
    #print(dictionary)

    corpus = [dictionary.doc2bow(text) for text in texts]
    corpora.MmCorpus.serialize(parameter.corpus, corpus)
    #print(corpus)

    print("----End step 2\n")
