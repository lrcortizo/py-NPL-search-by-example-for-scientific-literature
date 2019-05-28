import nltk

class Article:

    def __init__(self, pmid, article_title, abstract):
        self.pmid = pmid
        self.article_title = article_title
        self.abstract = abstract
        self.abstract_array = None

    def get_abstract_array(self):
        if self.abstract_array is None:
            punctuation = ".,?!:;(){}[]"

            self.abstract_array = nltk.sent_tokenize(self.abstract)
            
            for i, sentence in enumerate(self.abstract_array):
                for char in punctuation:
                    sentence = sentence.replace(char, ' %s ' % char)

                self.abstract_array[i] = sentence.split()

        return self.abstract_array
