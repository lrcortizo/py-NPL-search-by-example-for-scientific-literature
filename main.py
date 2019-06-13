import sys
import os
import getopt
import web_scraper
import text_processing
import text_similarities
from parameter import Parameter

def check_arguments(argv):
    #Initialize input params variables
    search_term = ''
    file = ''
    output_directory = ''
    max_results = ''

    #get params
    try:
        opts, args = getopt.getopt(argv,"hs:f:d:n:",["help","search=","file=","dir=","number="])
    except getopt.GetoptError:
        print ('main.py -s <search_term>  -f <input_file> -d <output_directory> -n <results_number>')
        sys.exit(2)

    #parse cli params
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print ('usage: main.py -s <search_term>  -f <input_file> -d <output_directory> -n <results_number>\n\
Options and arguments:\n\
  -s, --search: Search term to query pubmed\n\
  -f, --input_file: Input file to compare the results of query\n\
  -d, --dir: Output directory to the output and temporary files\n\
  -n, --results_number: Number of results in pubmed search')
            sys.exit()
        elif opt in ("-s", "--search"):
            search_term = arg
        elif opt in ("-f", "--file"):
            file = arg
            if not os.path.isfile(file) or not file.endswith('.txt'):
                print ('You must introduce a valid input file. Check path and format(txt)')
                sys.exit()
        elif opt in ("-d", "--dir"):
            output_directory = arg
        elif opt in ("-n", "--number"):
            max_results = arg

    #check valid params
    if search_term == '':
        print ('You must introduce a search term')
        sys.exit()
    elif file == '':
        print ('You must introduce a input file')
        sys.exit()
    elif output_directory == '':
        self.output_directory = './output/'
    elif max_results == '':
        self.max_results = '10'

    #object with the params configuration
    parameter = Parameter(search_term, file, output_directory, max_results)
    return parameter

if __name__ == "__main__":
    parameter = check_arguments(sys.argv[1:])
    #Step 1: Web scraping
    web_scraper.scrape(parameter)
    #Step 2: Text processing
    text_processing.process_text(parameter)
    #Step 3: Text similarities
    text_similarities.similarity(parameter)
