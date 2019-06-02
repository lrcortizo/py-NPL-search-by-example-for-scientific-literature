import os
import re

class Parameter:
    def __init__(self, search_term, file, output_directory, max_results):
        self.search_term = search_term
        self.file = file

        if output_directory is None:
            self.output_directory = './output/'
        else:
            self.output_directory = self.check_directory_format(output_directory)

        if max_results is None:
            self.max_results = '10'
        else:
            self.max_results = max_results

        self.scrapper_result = self.output_directory + "web_scrapper_results.xml"
        self.dictionary = self.output_directory + "dictionary.dict"
        self.corpus = self.output_directory + "corpus.mm"

    def check_directory_format(dir):
        pattern = re.compile(".+\/$")
        if pattern.match(dir) is None:
            dir  =  dir + "/"
        return dir

    def create_output_directory():
        if not os.path.exists(self.output_directory):
            os.mkdir(self.output_directory)
