import natural_language_processing

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
        self.abstract_array = natural_language_processing.preprocessing(self.abstract)
