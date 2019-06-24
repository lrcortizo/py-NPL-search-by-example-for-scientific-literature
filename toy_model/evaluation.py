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

def lsi_model(doc_clean,number_of_topics,words):
    dictionary,doc_term_matrix=prepare_corpus(doc_clean)
    tfidf = models.TfidfModel(doc_term_matrix)
    corpus_tfidf = tfidf[doc_term_matrix]
    lsi_model = models.LsiModel(corpus_tfidf, num_topics=number_of_topics, id2word = dictionary)

    for topic in lsi_model.print_topics(number_of_topics):
        print(topic)

    corpus_lsi = lsi_model[corpus_tfidf]

    for i, doc in enumerate(corpus_lsi):
        if i==0:
            print("\nSPORT\n")
        elif i==5:
            print("\nPOLITICS\n")
        elif i==10:
            print("\nBIOLOGY\n")
        elif i==15:
            print("\nMECHANICS\n")
        print(doc)

def compute_coherence_values(dictionary, doc_term_matrix, doc_clean, start, stop, step):
    coherence_values = []
    model_list = []
    for num_topics in range(start, stop, step):
        model = models.LsiModel(doc_term_matrix, num_topics=num_topics, id2word = dictionary)
        model_list.append(model)
        coherencemodel = models.CoherenceModel(model=model, texts=doc_clean, dictionary=dictionary, coherence='c_v')
        coherence_values.append(coherencemodel.get_coherence())
    return model_list, coherence_values

def plot_graph(texts_array, start, stop, step):
    dictionary, doc_term_matrix=prepare_corpus(texts_array)
    model_list, coherence_values = compute_coherence_values(dictionary, doc_term_matrix,texts_array,
                                                            start, stop, step)
    # Show graph
    x = range(start, stop, step)
    plt.plot(x, coherence_values)
    plt.xlabel("Number of topics")
    plt.ylabel("Coherence score")
    plt.legend(("coherence_values"), loc='best')
    plt.show()

words=10
number_of_topics=4
document_list=load_data("texts.txt")
texts_array=preprocess_data(document_list)
lsi_model(texts_array,number_of_topics,words)
#start,stop,step=2,20,1
#plot_graph(texts_array,start,stop,step)
