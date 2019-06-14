import gensim.similarities as similarities
from gensim import corpora
from gensim import models

def similarity(parameter):
    #load dictionary and corpus
    dictionary = corpora.Dictionary.load(parameter.dictionary)
    corpus_list = corpora.MmCorpus(parameter.corpus)
    corpus = list(corpus_list)

    tfidf = models.TfidfModel(list(corpus))
    corpus_tfidf = tfidf[corpus]

    lsi = models.LsiModel(corpus_tfidf, id2word=dictionary, num_topics=2)
    corpus_lsi = lsi[corpus_tfidf]

    vec_bow = dictionary.doc2bow(parameter.get_input_file_array())
    vec_lsi = lsi[vec_bow]

    index = similarities.MatrixSimilarity(corpus_lsi)
    sims = sorted(enumerate(index[vec_lsi]), key=lambda item: -item[1])
    print(sims)
