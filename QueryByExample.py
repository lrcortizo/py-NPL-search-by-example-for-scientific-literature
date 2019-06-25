import sys
import os
import re
import getopt
import logging
import multiprocessing
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
    if pattern.match(dir) is not None:
        return True
    return False

"""
Validate given results number/topics number
Input  : results number/topics number
Output : boolean check
"""
def valid_number(num):
    pattern = re.compile(r'\d+')
    if not pattern.match(num) is None and int(num) > 1:
        return True
    return False

######### parsing functions
# Help
def parse_help(parameter):
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
-c, --coherence_model: Coherence Model mode\n\
-h, --help: Shows the different options and arguments')
    sys.exit()

# Search term
def parse_search_term(arg, parameter):
    print(arg)
    if arg.strip():
        parameter.search_term = arg
    else:
        print ('You must introduce a valid search term. At least 3 characters.')
        sys.exit(1)

def parse_example_file(arg, parameter):
    if valid_file(arg):
        parameter.file = arg
    else:
        sys.exit(1)

# Output directory
def parse_output_directory(arg, parameter):
    if arg.strip() and valid_directory(arg):
        parameter.output_directory = arg
    else:
        print ('You must introduce a valid directory name.')
        sys.exit(1)

# Results number
def parse_max_results(arg, parameter):
    if valid_number(arg):
        parameter.max_results = arg
    else:
        print ('You must introduce a valid result number.')
        sys.exit(1)

# Coherence model mode
def parse_verbose(parameter):
    parameter.verbose = True

# Number of CPU processors
def parse_processors(arg, parameter):
    if int(arg) <= multiprocessing.cpu_count():
        parameter.processors = int(arg)
    else:
        print ('Too much processors. This computer has '+str(multiprocessing.cpu_count())+' processors')
        sys.exit(1)

# Number of topics
def parse_topics(arg, parameter):
    if valid_number(arg):
        parameter.topics = int(arg)
    else:
        print ('You must introduce a valid topics number.')
        sys.exit(1)

# Max topics
def parse_max_topics(arg, parameter):
    if valid_number(arg):
        parameter.max_topics = int(arg)
    else:
        print ('You must introduce a valid maximum topics number.')
        sys.exit(1)

# Coherence model mode
def parse_coherence_model(parameter):
    parameter.coherence_model = True

"""
Parse and validate input paramas
Input  : script params
Output : parameter object
"""
def cli_params(argv):
    """
    object with default params configuration
        search_term = '', file = '', output_directory = './',
        max_results = '10', verbose = False, processors = 1,
        topics =  20, max_topics = 200, coherence_model = False,
    """
    parameter = Parameter('', '', './', 10, False, 1, 20, 200, False)

    #dict to simulate switch case
    params = {'-h' : parse_help,
        '--help' : parse_help,
        '-s' : parse_search_term,
        '--search' : parse_search_term,
        '-f' : parse_example_file,
        '--file' : parse_example_file,
        '-d' : parse_output_directory,
        '--dir' : parse_output_directory,
        '-n' : parse_max_results,
        '--number' : parse_max_results,
        '-v' : parse_verbose,
        '--verbose' : parse_verbose,
        '-p' : parse_processors,
        '--procesors' : parse_processors,
        '-t' : parse_topics,
        '--topics' : parse_topics,
        '-m' : parse_max_topics,
        '--max_topics' : parse_max_topics,
        '-c' : parse_coherence_model,
        '--coherence_model'  : parse_coherence_model,
    }

    # get params
    try:
        opts, args = getopt.getopt(argv,"hs:f:d:n:vp:t:m:c",["help","search=","file=","dir=",
            "number=","verbose","processors=","topics=","max_topics=","coherence_model"])
    except getopt.GetoptError:
        print ('python QueryByExample.py -s <search_term> -f <input_file>')
        sys.exit(2)

    # parse and check cli params
    for opt, arg in opts:
        if arg:
            params[opt](arg, parameter)
        else:
            params[opt](parameter)

    # check required params
    if not parameter.search_term:
        print ('You must introduce a search term.')
        sys.exit(1)
    elif not parameter.file:
        print ('You must introduce a input file.')
        sys.exit(1)

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
    pmids, results = extract_data.extract(parameter)
    print_step("End step 1")

    # Step 2: Text processing
    print_step("Step 2: Natural language processing")
    articles, dicctionary, corpus = natural_language_processing.process_docs(parameter, pmids, results)
    print_step("End step 2")

    # Step 3: Text similarities
    print_step("Step 3: Text similarities")
    text_similarities.similarity(parameter, articles, dicctionary, corpus)
    print_step("End step 3")
