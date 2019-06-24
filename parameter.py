import os
import re
import string
from collections import defaultdict
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
from nltk.stem.porter import PorterStemmer

class Parameter:
    def __init__(self, search_term, file, output_directory, max_results,
        verbose, processors, topics, max_topics, coherence_model):
        self.search_term = search_term
        self.file = file
        self.output_directory = self.check_directory_format(output_directory)
        self.create_output_directory()
        self.max_results = max_results
        self.data_extraction_result = self.output_directory + "data_extraction_result.xml"
        self.final_result = self.output_directory + "similarities_result.txt"
        self.dictionary = self.output_directory + "dictionary.dict"
        self.corpus = self.output_directory + "corpus.mm"
        self.index = self.output_directory + "similarity.index"
        self.input_file_text = None
        self.verbose = verbose
        self.processors = processors
        self.topics = topics
        self.max_topics = max_topics
        self.coherence_model = coherence_model

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
            stemmer = PorterStemmer()
            #tokenize, remove stopwords and punctuation
            tokens = word_tokenize(text)
            tokens = [w.lower() for w in tokens]
            table = str.maketrans('', '', string.punctuation)
            stripped = [w.translate(table) for w in tokens]
            cleaned_tokens = [word for word in stripped if word.isalpha()]
            stop_words = set(stopwords.words('english'))
            stopped_tokens = [w for w in cleaned_tokens if not w in stop_words]
            stemmed_tokens = [stemmer.stem(w) for w in stopped_tokens]
            frequency = defaultdict(int)
            for word in stemmed_tokens:
                frequency[word] += 1
            self.input_file_text = [word for word in stemmed_tokens if frequency[word] > 1]


        return self.input_file_text
