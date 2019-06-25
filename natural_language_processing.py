import sys
import nltk
from bs4 import BeautifulSoup
from gensim import corpora
from article import Article

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
Parse xml results to a list of Article objects
Input  : path to xml, result PubMed IDs
Output : Article objects list
"""
def parse_xml(pmids, results):
    article_list = []
    # parse xml with BeautifoulSoup
    soup = BeautifulSoup(results,'xml')
    try:
        # parse xml file into a list of Article objects
        titles = soup.find_all('ArticleTitle')
        abstracts = soup.find_all('Abstract')
        print("* Parsing xml file...")
        for i in range(0, len(titles)):
            article_list.append(Article(pmids[i], titles[i].text, abstracts[i].text))
    except:
        print(20*"*"+" An error ocurred while parsing xml file. Try it again")
        sys.exit(1)
    return article_list

"""
Preprocessing of abstract texts before create corpus
Input  : Article list
Output : tokenized Article objects list
"""
def preprocessing(article_list):
    tokenized_list = []
    print("** Preprocessing and tokenizing texts...")
    for article in article_list:
        tokenized_list.append(article.get_abstract_array())

    return tokenized_list

"""
Creates dictionary and BOW corpus representation
Input  : parameter
         tokenized docs
Output : dicctionary
         corpus
"""
def prepare_corpus(parameter, tokenized_list):
    # Creates dictionary
    print("*** Building dictionary...")
    dictionary = corpora.Dictionary(tokenized_list)
    dictionary.save(parameter.dictionary)
    print ("**** Dictionary stored in " + parameter.dictionary)

    # Creating Bag of Words model
    print("***** Building BOW corpus...")
    corpus = [dictionary.doc2bow(array) for array in tokenized_list]
    corpora.MmCorpus.serialize(parameter.corpus, corpus)
    print ("****** Corpus stored in " + parameter.corpus)

    return dictionary, corpus

"""
Main method of natural_language_processing
Input  : parameter object, pmids list
Output : Article object list
         Tokenized text
         Dictionary
         Corpus
"""
def process_docs(parameter, pmids, results):
    # Check if necessary resources are avaliable
    nltk_check()

    # Parse xml file
    article_list = parse_xml(pmids, results)
    tokenized_list = preprocessing(article_list)

    # Creates dictionary and corpus
    dictionary, corpus = prepare_corpus(parameter, tokenized_list)

    return article_list, tokenized_list, dictionary, corpus
