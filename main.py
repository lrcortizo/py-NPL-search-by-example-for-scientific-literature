import sys
import os
import re
import getopt
import logging
import extract_data
import natural_language_processing
import text_similarities
from parameter import Parameter

def valid_file(file):
    if os.path.isfile(file) and file.endswith('.txt'):
        return True
    return False

def valid_directory(dir):
    pattern = re.compile("^[\w\_.@()-]+/?$")
    if pattern.match(dir) is None:
        return False
    return True

def valid_number(num):
    pattern = re.compile(r'\d+')
    if not pattern.match(num) is None and int(num) > 0:
        return True
    return False

def cli_params(argv):
    #Initialize input params variables
    search_term = ''
    file = ''
    output_directory = './' #default
    max_results = '10' #default
    verbose = False #default

    #get params
    try:
        opts, args = getopt.getopt(argv,"hs:f:d:n:v",["help","search=","file=","dir=","number=","verbose"])
    except getopt.GetoptError:
        print ('main.py -s <search_term>  -f <input_file> -d <output_directory> -n <results_number>')
        sys.exit(2)

    #parse and check cli params
    for opt, arg in opts:
        #Help
        if opt in ("-h", "--help"):
            print ('usage: main.py -s <search_term> -f <input_file> -d <output_directory> -n <results_number>\n\
Options and arguments:\n\
  -s, --search: Search term to query pubmed\n\
  -f, --input_file: Input file to compare the results of query\n\
  -d, --dir: Output directory to the output and temporary files\n\
  -n, --results_number: Number of results in pubmed search\n\
  -v, --verbose: Verbose mode')
            sys.exit()

        #Search term
        elif opt in ("-s", "--search"):
            if len(arg) > 2:
                search_term = arg
            else:
                print ('You must introduce a valid search term. At least 3 characters.')
                sys.exit()

        #Reference file
        elif opt in ("-f", "--file"):
            if valid_file(arg):
                file = arg
            else:
                print ('You must introduce a valid input file. Check path and format(txt).')
                sys.exit()

        #Output directory
        elif opt in ("-d", "--dir"):
            if arg.strip() and valid_directory(arg):
                output_directory = arg
            else:
                print ('You must introduce a valid directory name.')
                sys.exit()

        #Results number
        elif opt in ("-n", "--number"):
            if valid_number(arg):
                max_results = arg
            else:
                print ('You must introduce a valid result number.')
                sys.exit()

        #Verbose mode
        elif opt in ("-v", "--verbose"):
            verbose = True

    #check required params
    if not search_term:
        print ('You must introduce a search term.')
        sys.exit()
    elif not file:
        print ('You must introduce a input file.')
        sys.exit()

    #object with the params configuration
    parameter = Parameter(search_term, file, output_directory, max_results, verbose)
    return parameter

def print_step(msg):
    num = int((100-(len(msg)+2))/2)
    print("\n"+100*"*"+"\n"+num*"-"+" "+msg+" "+num*"-"+"\n"+100*"*"+"\n")


if __name__ == "__main__":
    parameter = cli_params(sys.argv[1:])

    if parameter.verbose:
        logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

    #Step 1: Web scraping
    print_step("Step 1: Scraping pubmed database")
    extract_data.extract(parameter)
    print_step("End step 1")

    #Step 2: Text processing
    print_step("Step 2: Natural language processing")
    articles = natural_language_processing.process_docs(parameter)
    print_step("End step 2")

    #Step 3: Text similarities
    print_step("Step 3: Text similarities")
    text_similarities.similarity(parameter, articles)
    print_step("End step 3")
