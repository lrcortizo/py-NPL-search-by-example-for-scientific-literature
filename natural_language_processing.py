import sys
import string
import nltk
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
from nltk.stem.porter import PorterStemmer
from bs4 import BeautifulSoup
from gensim import corpora
from collections import defaultdict
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

def frequency(words):
    frequency = defaultdict(int)
    for word in words:
        frequency[word] += 1

    words = [word for word in words if frequency[word] > 1]
    return words

"""
Preprocessing of texts
Input  : plain text
Output : tokenized text
"""
def preprocessing(text):
    tokenized_list = []
    stemmer = PorterStemmer()
    table = str.maketrans('', '', string.punctuation)
    stop_words = set(stopwords.words('english'))

    # split into words
    tokens = word_tokenize(text)
    # convert to lower case
    tokens = [w.lower() for w in tokens]
    # remove punctuation from each word
    stripped = [w.translate(table) for w in tokens]
    # remove remaining tokens that are not alphabetic
    cleaned_tokens = [word for word in stripped if word.isalpha()]
    # filter out stop words
    stopped_tokens = [w for w in cleaned_tokens if not w in stop_words]
    # stemming
    stemmed_tokens = [stemmer.stem(w) for w in stopped_tokens]
    # remove frequency 1 words
    tokenized_text = frequency(stemmed_tokens)

    return tokenized_text

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
         Dictionary
         Corpus
"""
def process_docs(parameter, pmids, results):
    # Check if necessary resources are avaliable
    nltk_check()

    # Parse xml file
    article_list = parse_xml(pmids, results)

    # Tokenize example file
    file = open(parameter.file, mode='r')
    text = file.read()
    file.close()
    parameter.input_file_text = preprocessing(text)

    # tokenized article texts
    tokenized_list = [article.abstract_array for article in article_list]

    # Creates dictionary and corpus
    dictionary, corpus = prepare_corpus(parameter, tokenized_list)

    return article_list, dictionary, corpus
