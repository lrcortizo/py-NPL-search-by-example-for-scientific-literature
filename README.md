## Query-by-example-for-systematic-scientific-literature-searches
##### Query-by-example pipeline for scientific literature searches using web scraping, natural language processing and machine learning clasifiers
###### Usage:
```
main.py -s <search_term>  -f <input_file> -d <output_directory> -n <results_number>
Options and arguments:
  -s, --search: Search term to query pubmed
  -f, --input_file: Input file to compare the results of query
  -d, --dir: Output directory to the output and temporary files
  -n, --results_number: Number of results in pubmed search
```
###### Requirements:
- Python >= 3.7.2
  - biopython
  - gensim
  - BeautifulSoup
  - nltk
  - os
  - re
  - sys
  - getopt
- C Compiler

