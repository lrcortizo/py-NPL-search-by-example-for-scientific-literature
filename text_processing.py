import gensim
import gensim.similarities as similarities
import logging
from bs4 import BeautifulSoup
from article import Article
from gensim import corpora
from gensim import models

def parse_xml(xml_path):
    try:
        article_list = []

        infile = open(xml_path,"r")
        contents = infile.read()
        soup = BeautifulSoup(contents,'xml')

        pmids = soup.find_all('PMID')
        titles = soup.find_all('ArticleTitle')
        abstracts = soup.find_all('Abstract')

        for i in range(0, len(pmids)):
            article_list.append(Article(pmids[i].text, titles[i].text, abstracts[i].text))
        return article_list
    except:
        return None

if __name__ == '__main__':
    #logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

    xml_path = 'tmp/web_scrapper_results.xml'
    article_list = parse_xml(xml_path)
    #print(article_list[0].get_abstract_array())
    texts = []
    for article in article_list:
        texts.append(article.get_abstract_array())
    #print(texts)

    dictionary = corpora.Dictionary(texts)
    dictionary.save('tmp/dictionary.dict')
    #print(dictionary)

    corpus = [dictionary.doc2bow(text) for text in texts]
    corpora.MmCorpus.serialize('tmp/corpus.mm', corpus)
    #print(corpus)

    tfidf = models.TfidfModel(corpus)
    corpus_tfidf = tfidf[corpus]
    #for doc in corpus_tfidf:
    #    print(doc)
    lsi = models.LsiModel(corpus_tfidf, id2word=dictionary, num_topics=2)
    corpus_lsi = lsi[corpus_tfidf]
    #for doc in corpus_lsi:
        #print(doc)

    lsi = models.LsiModel(corpus, id2word=dictionary, num_topics=2)
    vec_bow = dictionary.doc2bow(article_list[0].abstract.lower().split())
    vec_lsi = lsi[vec_bow]


    index = similarities.MatrixSimilarity(lsi[corpus])
    sims = sorted(enumerate(index[vec_lsi]), key=lambda item: -item[1])
    #print(sims)
