import gensim.similarities as similarities
from gensim import corpora
from gensim import models

def build_output(sims, articles):
    results = "\n"+20*"-"+"RESULTS"+20*"-"+"\n\n"
    for x in sims:
        results += articles[x[0]].article_title
        results += "\nAccuracy: "+str(x[1])+"\n\n"

    return results

def write_output(parameter, results):
    #write results into a txt file
    f=open(parameter.output_result ,"w", encoding='utf-8')
    f.write(results)
    f.close()
    print ("Results stored in " + parameter.output_result)


def similarity(parameter, articles):
    #load dictionary and corpus
    dictionary = corpora.Dictionary.load(parameter.dictionary)
    corpus_list = corpora.MmCorpus(parameter.corpus)
    corpus = list(corpus_list)

    #Creating tfidf model
    tfidf = models.TfidfModel(corpus)
    corpus_tfidf = tfidf[corpus]
    #Double wrapping with lsi
    lsi = models.LsiModel(corpus_tfidf, id2word=dictionary, num_topics=20)
    lsi.print_topics(20)
    corpus_lsi = lsi[corpus_tfidf]

    #Reference file to compare
    vec_bow = dictionary.doc2bow(parameter.get_input_file_array())
    vec_tfidf = tfidf[vec_bow]
    vec_lsi = lsi[vec_tfidf]

    #Similarites
    index = similarities.MatrixSimilarity(corpus_lsi)
    sims = sorted(enumerate(index[vec_lsi]), key=lambda item: -item[1])
    #Print sorted articles
    results = build_output(sims, articles)
    print(results)
    write_output(parameter, results)
