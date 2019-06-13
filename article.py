import string
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords

class Article:

    def __init__(self, pmid, article_title, abstract):
        self.pmid = pmid
        self.article_title = article_title
        self.abstract = abstract
        self.abstract_array = None

    def get_abstract_array(self):
        #tokenize the abstract
        if self.abstract_array is None:
            # split into words
            tokens = word_tokenize(self.abstract)
            # convert to lower case
            tokens = [w.lower() for w in tokens]
            # remove punctuation from each word
            table = str.maketrans('', '', string.punctuation)
            stripped = [w.translate(table) for w in tokens]
            # remove remaining tokens that are not alphabetic
            words = [word for word in stripped if word.isalpha()]
            # filter out stop words
            stop_words = set(stopwords.words('english'))
            self.abstract_array = [w for w in words if not w in stop_words]

        return self.abstract_array
