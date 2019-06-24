import sys
import os
import re
import getopt
import logging
import extract_data
import natural_language_processing
import text_similarities
from parameter import Parameter

"""
Validate input example file
Input  : input file
Output : boolean check
"""
def valid_file(file):
    if os.path.isfile(file):
        if os.access(file, os.R_OK):
            if file.endswith('.txt'):
                if os.path.getsize(file) > 0:
                    return True
                else:
                    print ('The example file is empty.')
            else:
                print ('The example file must have txt format.')
        else:
            print ('The example file must have read permissions.')
    else:
        print ('Example file does not exist.')
    return False

"""
Validate given output directory
Input  : output directory path
Output : boolean check
"""
def valid_directory(dir):
    pattern = re.compile("^[\w\_.@()-]+/?$")
    if pattern.match(dir) is None:
        return False
    return True

"""
Validate given results number
Input  : results number
Output : boolean check
"""
def valid_number(num):
    pattern = re.compile(r'\d+')
    if not pattern.match(num) is None and int(num) > 0:
        return True
    return False

"""
Parse and validate input paramas
Input  : script params
Output : parameter object
"""
def cli_params(argv):
    # Initialize input params variables
    search_term = ''
    file = ''
    output_directory = './' # default
    max_results = '10' # default
    verbose = False # default
    processors = 1 # default
    topics =  20 # default
    max_topics = 200 # default
    coherence_model = False # default

    # get params
    try:
        opts, args = getopt.getopt(argv,"hs:f:d:n:vp:t:m:c",["help","search=","file=","dir=",
            "number=","verbose","processors","topics","max_topics","coherence_model"])
    except getopt.GetoptError:
        print ('python QueryByExample.py -s <search_term> -f <input_file> ')
        sys.exit(2)

    # parse and check cli params
    for opt, arg in opts:
        # Help
        if opt in ("-h", "--help"):
            print ('usage: python QueryByExample.py -s <search_term> -f <input_file> \n\
Options and arguments:\n\
  -s, --search: Search term to query pubmed\n\
  -f, --input_file: Input file to compare the results of query\n\
  -d, --dir: Output directory to the output and temporary files\n\
  -n, --results_number: Number of results in pubmed search\n\
  -v, --verbose: Verbose mode\n\
  -p, --procesors: Number of CPU processors\n\
  -t, --topics: Number of topics\n\
  -m, --max_topics: Maximum number of topics\n\
  -c, --coherence_model: Coherence Model mode')
            sys.exit()

        # Search term
        elif opt in ("-s", "--search"):
            if arg.strip():
                search_term = arg
            else:
                print ('You must introduce a valid search term. At least 3 characters.')
                sys.exit()

        # Example file
        elif opt in ("-f", "--file"):
            if valid_file(arg):
                file = arg
            else:
                sys.exit()

        # Output directory
        elif opt in ("-d", "--dir"):
            if arg.strip() and valid_directory(arg):
                output_directory = arg
            else:
                print ('You must introduce a valid directory name.')
                sys.exit()

        # Results number
        elif opt in ("-n", "--number"):
            if valid_number(arg):
                max_results = arg
            else:
                print ('You must introduce a valid result number.')
                sys.exit()

        # Number of CPU processors
        elif opt in ("-p", "--procesors"):
            processors = arg

        # Number of topics
        elif opt in ("-t", "--topics"):
            topics = arg

        # Max topics
        elif opt in ("-m", "--max_topics"):
            max_topics = arg

        # Coherence model mode
        elif opt in ("-c", "--coherence_model"):
            coherence_model = True

    # check required params
    if not search_term:
        print ('You must introduce a search term.')
        sys.exit()
    elif not file:
        print ('You must introduce a input file.')
        sys.exit()

    # object with the params configuration
    parameter = Parameter(search_term, file, output_directory, max_results, verbose,
        processors, topics, max_topics, coherence_model)
    return parameter

"""
Prints info messages when start and finish each step
Input  : message to print
"""
def print_step(msg):
    num = int((100-(len(msg)+2))/2)
    print("\n"+100*"*"+"\n"+num*"-"+" "+msg+" "+num*"-"+"\n"+100*"*"+"\n")

"""
Main class of the application
Input  : input script parameters
"""
if __name__ == "__main__":
    parameter = cli_params(sys.argv[1:])
    # full logging if verbose mode is on
    if parameter.verbose:
        logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

    # Step 1: Web scraping
    print_step("Step 1: Scraping pubmed database")
    extract_data.extract(parameter)
    print_step("End step 1")

    # Step 2: Text processing
    print_step("Step 2: Natural language processing")
    articles = natural_language_processing.process_docs(parameter)
    print_step("End step 2")

    # Step 3: Text similarities
    print_step("Step 3: Text similarities")
    text_similarities.similarity(parameter, articles)
    print_step("End step 3")
