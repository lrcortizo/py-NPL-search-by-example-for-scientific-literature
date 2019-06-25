## Query-by-example-for-systematic-scientific-literature-searches
##### Query-by-example tool for scientific literature searches using web scraping, natural language processing and machine learning
###### Usage:
```
python QueryByExample.py -s <search_term>  -f <input_file>
Options and arguments:
  -s, --search: Search term to query pubmed
  -f, --input_file: Input file to compare the results of query
  -d, --dir: Output directory to the output and temporary files (Optional, current directory by default)
  -n, --results_number: Number of results in pubmed search (Optional, 10 results by default)
  -v, --verbose: Verbose mode (Optional, no verbose mode by default)
  -p, --processors: Number of CPU processors (Optional, 1 processor by default)
  -t, --topics: Number of topics (Optional, 20 topics by default)
  -m, --max_topics: Maximum number of topics (Optional, 200 by default)
  -c, --coherence_model: Coherence Model mode (Optional, no coherence model mode by default)
  -h, --help: Shows the different options and arguments
```
###### Requirements:
- Python >= 3.7.2
  - biopython: ```pip install biopython```
  - gensim: ```pip install gensim```
  - BeautifulSoup: ```pip install beautifulsoup4```
  - lxml: ```pip install lxml```
  - nltk: ```pip install nltk```

###### Notes:
- Scripts need execute permissions, as well as the input file requires read permissions and the output directory read/write permissions.
- The search term must be enclosed in quotation marks (if it consists of more than one word).
- The output directory cannot contain special characters, only _ - ( ) @ are allowed.

If some argument is incorrect, an error message will be printed in the console indicating which data is wrong. The script execution will generate the following files in the output directory:
- data_extraction_result.xml
- similarities_result.txt
- tmp/ directory
  - dictionary.dict
  - corpus.mm
  - similarity.index

The console output will show different messages depending on whether the verbose mode is activated. The output is divided in three 3 steps (data extraction, natural language processing and classification) showing the info of each one.
