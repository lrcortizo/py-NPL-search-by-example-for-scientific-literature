import nltk

class Article:

    def __init__(self, pmid, article_title, abstract):
        self.pmid = pmid
        self.article_title = article_title
        self.abstract = abstract
        self.abstract_array = None

    def get_abstract_array(self):
        #tokenize the abstract
        if self.abstract_array is None:
            stoplist = set('for a of the and to in'.split())
            self.abstract_array = [word for word in self.abstract.lower().split() if word not in stoplist]

        return self.abstract_array
