import sys
import os
import getopt
import web_scraper
import natural_language_processing
import text_similarities
from parameter import Parameter

def cli_params(argv):
    #Initialize input params variables
    search_term = ''
    file = ''
    output_directory = ''
    max_results = ''
    verbose = False

    #get params
    try:
        opts, args = getopt.getopt(argv,"hs:f:d:n:v",["help","search=","file=","dir=","number=","verbose"])
    except getopt.GetoptError:
        print ('main.py -s <search_term>  -f <input_file> -d <output_directory> -n <results_number>')
        sys.exit(2)

    #parse cli params
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print ('usage: main.py -s <search_term> -f <input_file> -d <output_directory> -n <results_number>\n\
Options and arguments:\n\
  -s, --search: Search term to query pubmed\n\
  -f, --input_file: Input file to compare the results of query\n\
  -d, --dir: Output directory to the output and temporary files\n\
  -n, --results_number: Number of results in pubmed search\n\
  -v, --verbose: Verbose mode')
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
        elif opt in ("-v", "--verbose"):
            verbose = True

    #check valid params
    if not search_term.strip():
        print ('You must introduce a search term')
        sys.exit()
    elif not file.strip():
        print ('You must introduce a input file')
        sys.exit()
    elif not output_directory.strip():
        self.output_directory = './'
    elif not max_results.strip():
        self.max_results = '10'

    #object with the params configuration
    parameter = Parameter(search_term, file, output_directory, max_results, verbose)
    return parameter

if __name__ == "__main__":
    parameter = cli_params(sys.argv[1:])
    #Step 1: Web scraping
    web_scraper.scrape(parameter)
    #Step 2: Text processing
    articles = natural_language_processing.process_text(parameter)
    #Step 3: Text similarities
    text_similarities.similarity(parameter)
