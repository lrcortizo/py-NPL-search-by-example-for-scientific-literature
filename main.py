import sys
import getopt
import web_scrapper
import text_processing
from parameter import Parameter

def check_arguments(argv):
    search_term = ''
    file = ''
    output_directory = ''
    max_results = ''

    try:
        opts, args = getopt.getopt(argv,"hs:f:d:n:",["help=","search=","file=","dir=","number="])
    except getopt.GetoptError:
        print ('main.py -s <search_term>  -f <input_file> -d <output_directory> -n <results_number>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print ('main.py -s <search_term>  -f <input_file> -d <output_directory> -n <results_number>')
            sys.exit()
        elif opt in ("-s", "--search"):
            search_term = arg
        elif opt in ("-f", "--file"):
            file = arg
        elif opt in ("-d", "--dir"):
            output_directory = arg
        elif opt in ("-n", "--number"):
            max_results = arg

    parameter = Parameter(search_term, file, output_directory, max_results)
    return parameter

if __name__ == "__main__":
    parameter = check_arguments(sys.argv[1:])
