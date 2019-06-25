import string
from collections import defaultdict
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
from nltk.stem.porter import PorterStemmer

"""
Class that stores articles information
"""
class Article:
    """
    Constructor of Article class
    Input  :  PubMed ID
              article title
              abstract
    """
    def __init__(self, pmid, article_title, abstract):
        self.pmid = pmid
        self.article_title = article_title
        self.abstract = abstract
        self.abstract_array = None

    def frequency(self, words):
        frequency = defaultdict(int)
        for word in words:
            frequency[word] += 1

        words = [word for word in words if frequency[word] > 1]
        return words

    def get_abstract_array(self):
        #tokenize the abstract
        if self.abstract_array is None:
            stemmer = PorterStemmer()
            # split into words
            tokens = word_tokenize(self.abstract)
            # convert to lower case
            tokens = [w.lower() for w in tokens]
            # remove punctuation from each word
            table = str.maketrans('', '', string.punctuation)
            stripped = [w.translate(table) for w in tokens]
            # remove remaining tokens that are not alphabetic
            cleaned_tokens = [word for word in stripped if word.isalpha()]
            # filter out stop words
            stop_words = set(stopwords.words('english'))
            stopped_tokens = [w for w in cleaned_tokens if not w in stop_words]
            stemmed_tokens = [stemmer.stem(w) for w in stopped_tokens]

            self.abstract_array = self.frequency(stemmed_tokens)

        return self.abstract_array
