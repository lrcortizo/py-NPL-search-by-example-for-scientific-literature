import os
import re
import string
from collections import defaultdict
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
from nltk.stem.porter import PorterStemmer

"""
Class that stores the application configuration
"""
class Parameter:
    """
    Constructor of Paramter class
    Input  :  search term
              example file
              output directory
              results number
              verbose mode
              processors number
              number of topics
              maximum topics number
              choerence model mode
    """
    def __init__(self, search_term, file, output_directory, max_results,
        verbose, processors, topics, max_topics, coherence_model):
        self.search_term = search_term
        self.file = file
        self.output_directory = self.create_output_directory(output_directory)
        self.max_results = max_results
        self.data_extraction_result = self.output_directory + "data_extraction_result.xml"
        self.final_result = self.output_directory + "similarities_result.txt"
        self.dictionary = self.output_directory + "tmp/dictionary.dict"
        self.corpus = self.output_directory + "tmp/corpus.mm"
        self.index = self.output_directory + "tmp/similarity.index"
        self.input_file_text = None
        self.verbose = verbose
        self.processors = processors
        self.topics = topics
        self.max_topics = max_topics
        self.coherence_model = coherence_model

    """
    Checks format of output directory name
    Input  : script params
    Output : parameter object
    """
    def check_directory_format(self, dir):
        # Check format of directory path
        pattern = re.compile(".+\/$")
        if pattern.match(dir) is None:
            dir  =  dir + "/"
        return dir

    """
    Checks and creates output directory if it not exists
    Input  : script params
    Output : parameter object
    """
    def create_output_directory(self, output_directory):
        dir = self.check_directory_format(output_directory)
        try:
            # Create output directory if not exists
            if not os.path.exists(dir):
                all_dirs = dir + "tmp/"
                os.mkdirs(all_dirs)
            else:
                if os.access(dir, os.W_OK):
                    if not os.path.exists(dir+"tmp/"):
                        os.mkdir(dir+"tmp/")
                else:
                    print("The output directory does not have wrtie permissions")
                    sys.exit(1)
        except:
            print("Cannot create output directory, check write permissions")
            sys.exit(1)
        return dir

    def get_input_file_array(self):
        if self.input_file_text is None:
            # parse file
            file = open(self.file, mode='r')
            text = file.read()
            file.close()
            stemmer = PorterStemmer()
            # tokenize, remove stopwords and punctuation
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
