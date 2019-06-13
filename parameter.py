import os
import re
import string
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords

class Parameter:
    def __init__(self, search_term, file, output_directory, max_results):
        self.search_term = search_term
        self.file = file
        self.output_directory = self.check_directory_format(output_directory)
        self.max_results = max_results
        self.scrapper_result = self.output_directory + "web_scraper_results.xml"
        self.dictionary = self.output_directory + "dictionary.dict"
        self.corpus = self.output_directory + "corpus.mm"
        self.input_file_text = None

    def check_directory_format(self, dir):
        #Check format of directory path
        pattern = re.compile(".+\/$")
        if pattern.match(dir) is None:
            dir  =  dir + "/"
        return dir

    def create_output_directory(self):
        #Create output directory if not exists
        if not os.path.exists(self.output_directory):
            os.mkdir(self.output_directory)

    def get_input_file_array(self):
        if self.input_file_text is None:
            #parse file
            file = open(self.file, mode='r')
            text = file.read()
            file.close()

            #tokenize, remove stopwords and punctuation
            tokens = word_tokenize(text)
            tokens = [w.lower() for w in tokens]
            table = str.maketrans('', '', string.punctuation)
            stripped = [w.translate(table) for w in tokens]
            words = [word for word in stripped if word.isalpha()]
            stop_words = set(stopwords.words('english'))
            self.input_file_text = [w for w in words if not w in stop_words]

        return self.input_file_text
