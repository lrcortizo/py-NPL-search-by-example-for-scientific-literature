import gensim.similarities as similarities
from gensim import corpora
from gensim import models

def similarity(parameter):
    print("----Step 3: Text similarities")
    #load dictionary and corpus
    dictionary = corpora.Dictionary.load(parameter.dictionary)
    corpus_list = corpora.MmCorpus(parameter.corpus)
    corpus = list(corpus_list)

    tfidf = models.TfidfModel(list(corpus))
    corpus_tfidf = tfidf[corpus]
    #for doc in corpus_tfidf:
    #    print(doc)
    lsi = models.LsiModel(corpus_tfidf, id2word=dictionary, num_topics=2)
    corpus_lsi = lsi[corpus_tfidf]
    #for doc in corpus_lsi:
        #print(doc)

    vec_bow = dictionary.doc2bow(parameter.get_input_file_array())
    vec_lsi = lsi[vec_bow]

    index = similarities.MatrixSimilarity(corpus_lsi)
    sims = sorted(enumerate(index[vec_lsi]), key=lambda item: -item[1])
    print(sims)

    print("----End step 3\n")
