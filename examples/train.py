import os.path
import string
from gensim import corpora
from gensim import models
from collections import defaultdict
from nltk.tokenize import word_tokenize
from nltk.stem.porter import PorterStemmer
from nltk.corpus import stopwords
import matplotlib.pyplot as plt

def load_data(file_name):
    documents_list = []
    with open(file_name ,"r") as f:
        for line in f.readlines():
            text = line.strip()
            documents_list.append(text)
    print("Number of docs: "+str(len(documents_list)))
    return documents_list

def preprocess_data(doc_set):
    stemmer = PorterStemmer()
    stop_words = set(stopwords.words('english'))
    table = str.maketrans('', '', string.punctuation)
    texts = []
    for i in doc_set:
        tokens = word_tokenize(i)
        tokens_lower = [w.lower() for w in tokens]
        stripped = [w.translate(table) for w in tokens_lower]
        cleaned_tokens = [word for word in stripped if word.isalpha()]
        stopped_tokens = [w for w in cleaned_tokens if not w in stop_words]
        stemmed_tokens = [stemmer.stem(w) for w in stopped_tokens]
        frequency = defaultdict(int)
        for word in stemmed_tokens:
            frequency[word] += 1
        final_tokens = [word for word in stemmed_tokens if frequency[word] > 1]
        texts.append(final_tokens)
    return texts

def prepare_corpus(doc_clean):
    dictionary = corpora.Dictionary(doc_clean)
    doc_term_matrix = [dictionary.doc2bow(doc) for doc in doc_clean]
    return dictionary,doc_term_matrix

def create_gensim_lsa_model(doc_clean,number_of_topics,words):

    dictionary,doc_term_matrix=prepare_corpus(doc_clean)
    tfidf = models.TfidfModel(doc_term_matrix)
    corpus_tfidf = tfidf[doc_term_matrix]
    # generate LSA model
    lsamodel = models.LsiModel(corpus_tfidf, num_topics=number_of_topics, id2word = dictionary)  # train model
    print(lsamodel.print_topics(number_of_topics))
    corpus_lsa = lsamodel[corpus_tfidf]

    for doc in corpus_lsa:
        print(doc)

    return lsamodel

def compute_coherence_values(dictionary, doc_term_matrix, doc_clean, stop, start=2, step=3):
    coherence_values = []
    model_list = []
    for num_topics in range(start, stop, step):
        # generate LSA model
        model = models.LsiModel(doc_term_matrix, num_topics=num_topics, id2word = dictionary)  # train model
        model_list.append(model)
        coherencemodel = models.CoherenceModel(model=model, texts=doc_clean, dictionary=dictionary, coherence='c_v')
        coherence_values.append(coherencemodel.get_coherence())
    return model_list, coherence_values

def plot_graph(doc_clean,start, stop, step):
    dictionary,doc_term_matrix=prepare_corpus(doc_clean)
    model_list, coherence_values = compute_coherence_values(dictionary, doc_term_matrix,doc_clean,
                                                            stop, start, step)
    # Show graph
    x = range(start, stop, step)
    plt.plot(x, coherence_values)
    plt.xlabel("Number of Topics")
    plt.ylabel("Coherence score")
    plt.legend(("coherence_values"), loc='best')
    plt.show()

# LSA Model
number_of_topics=4
words=10
document_list=load_data("texts.txt")
clean_text=preprocess_data(document_list)
model=create_gensim_lsa_model(clean_text,number_of_topics,words)
#start,stop,step=2,12,1
#plot_graph(clean_text,start,stop,step)
