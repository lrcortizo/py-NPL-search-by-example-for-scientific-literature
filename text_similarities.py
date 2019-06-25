import gensim.similarities as similarities
import natural_language_processing
from gensim import corpora
from gensim import models

"""
Compute coherence values to get the best model based on topics number
Input  : dicctionary
         document-term MatrixSimilarity
         tokenized texts
         start number of topics
         maximum number of topics
         increase step
Output : best scored model
         best coherence value
         best number of topics
"""
def get_coherence_value(dictionary, doc_term_matrix, tokenized_list, start, stop, step, processors):
    coherence_value = 0
    best_model = None
    best_num_topics = 0
    for num_topics in range(start, stop, step):
        # generate LSA model
        model = models.LsiModel(doc_term_matrix, num_topics=num_topics, id2word = dictionary)  # train model
        coherencemodel = models.CoherenceModel(model=model, texts=tokenized_list, dictionary=dictionary,
                                coherence='c_v', processes=processors)
        if(coherence_value < coherencemodel.get_coherence()):
            coherence_value = coherencemodel.get_coherence()
            best_model = model
            best_num_topics = num_topics
    return best_model, coherence_value, best_num_topics

"""
Builds the application output
Input  : similarities
         articles list
         increase step
Output : application output
"""
def build_output(sims, articles):
    results = "\n"+20*"-"+"RESULTS"+20*"-"+"\n\n"
    for x in sims:
        results += articles[x[0]].article_title
        results += "\nAccuracy: "+str(x[1])+"\n\n"

    return results

"""
Writes the application output in a text file
Input  : parameter object
         application output
"""
def write_output(parameter, results):
    #write results into a txt file
    f=open(parameter.final_result ,"w", encoding='utf-8')
    f.write(results)
    f.close()
    print ("Results stored in " + parameter.final_result)

"""
Main method of text_similarities
Input  : parameter object
         articles list
         tokenized texts
         dictionary
         corpus
"""
def similarity(parameter, articles, tokenized_list, dictionary, corpus):
    # Creating tfidf model
    tfidf = models.TfidfModel(corpus)
    doc_term_matrix = tfidf[corpus]

    if parameter.coherence_model:
        # Get best LSI model
        lsi_model, coherence_value, num_topics = get_coherence_value(dictionary, doc_term_matrix, tokenized_list,
                                                    2, parameter.max_topics, 1, parameter.processors)
    else:
        lsi_model = models.LsiModel(doc_term_matrix, num_topics=parameter.topics, id2word = dictionary)

    # Double wrapping with lsi
    corpus_lsi = lsi_model[doc_term_matrix]

    # Reference file to compare
    vec_bow = dictionary.doc2bow(parameter.get_input_file_array())
    vec_tfidf = tfidf[vec_bow]
    vec_lsi = lsi_model[vec_tfidf]

    # Similarites
    index = similarities.MatrixSimilarity(corpus_lsi)
    index.save(parameter.index)
    sims = sorted(enumerate(index[vec_lsi]), key=lambda item: -item[1])

    # Print sorted articles
    results = build_output(sims, articles)
    if parameter.coherence_model:
        print("\nCoherence value: "+str(coherence_value)+", Number of topics: "+str(num_topics))
    print(results)
    write_output(parameter, results)
