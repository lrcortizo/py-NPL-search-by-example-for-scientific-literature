import gensim.similarities as similarities
import natural_language_processing
from gensim import corpora
from gensim import models

def get_coherence_value(dictionary, doc_term_matrix, doc_clean, start, stop, step):
    coherence_value = 0
    best_model = None
    best_num_topics = 0
    for num_topics in range(start, stop, step):
        # generate LSA model
        model = models.LsiModel(doc_term_matrix, num_topics=num_topics, id2word = dictionary)  # train model
        coherencemodel = models.CoherenceModel(model=model, texts=doc_clean, dictionary=dictionary, coherence='c_v', processes=2)
        if(coherence_value < coherencemodel.get_coherence()):
            coherence_value = coherencemodel.get_coherence()
            best_model = model
            best_num_topics = num_topics
    return best_model, coherence_value, best_num_topics

def build_output(sims, articles):
    results = "\n"+20*"-"+"RESULTS"+20*"-"+"\n\n"
    for x in sims:
        results += articles[x[0]].article_title
        results += "\nAccuracy: "+str(x[1])+"\n\n"

    return results

def write_output(parameter, results):
    #write results into a txt file
    f=open(parameter.final_result ,"w", encoding='utf-8')
    f.write(results)
    f.close()
    print ("Results stored in " + parameter.final_result)

def similarity(parameter, articles):
    coherence_value = 0
    num_topics = 0
    #load dictionary and corpus
    doc_arrays = natural_language_processing.preprocessing(articles)
    dictionary = corpora.Dictionary.load(parameter.dictionary)
    corpus_list = corpora.MmCorpus(parameter.corpus)
    doc_term_matrix = list(corpus_list)

    #Creating tfidf model
    tfidf = models.TfidfModel(doc_term_matrix)
    corpus_tfidf = tfidf[doc_term_matrix]
    #Double wrapping with lsi

    #lsi_model, coherence_value, num_topics = get_coherence_value(dictionary, corpus_tfidf, doc_arrays, 10, 50, 1)
    lsi_model = models.LsiModel(corpus_tfidf, num_topics=50, id2word = dictionary)
    corpus_lsi = lsi_model[corpus_tfidf]

    #Reference file to compare
    vec_bow = dictionary.doc2bow(parameter.get_input_file_array())
    vec_tfidf = tfidf[vec_bow]
    vec_lsi = lsi_model[vec_tfidf]

    #Similarites
    index = similarities.MatrixSimilarity(corpus_lsi)
    index.save(parameter.index)
    index = similarities.MatrixSimilarity.load(parameter.index)
    sims = sorted(enumerate(index[vec_lsi]), key=lambda item: -item[1])
    #Print sorted articles
    results = build_output(sims, articles)
    if coherence_value > 0:
        print("\nCoherence value: "+str(coherence_value)+", Number of topics: "+str(num_topics))
    print(results)
    write_output(parameter, results)
