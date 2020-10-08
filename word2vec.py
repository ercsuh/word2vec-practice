from gensim.models import Word2Vec, KeyedVectors
import pandas as pd

class Sentences(object):
    """
    Extracts abstracts from Pubtator data.
    """
    def __init__(self, filename):
        self.filename = filename
 
    def __iter__(self):
        for line in open(self.filename):
            if "|a|" in line:
                yield line.lower().split("|")[2].split()
                # yield line

def create_word2vec(sentences):
    model = Word2Vec(sentences, size=200, window=5, min_count=1, workers=4)
    model.save("word2vec.model")
    return model

def get_gene_disease_pairs(filename):
    df = pd.read_csv(filename, sep="\t")
    pairs = df[["doid_name", "gene_symbol"]].values.tolist()
    return pairs

def get_scores(model, pairs):
    outfile = open("similarity_scores.tsv", "w")
    outfile.write("disease\tgene\tscore\n")
    for pair in pairs:
        if all(word.lower().split()[0] in model.wv.vocab for word in pair):
            print(f"pair exists in vocabulary: {pair[0].lower().split()[0]} AND {pair[1].lower().split()[0]}")
            score = model.wv.similarity(pair[0].lower().split()[0], pair[1].lower().split()[0])
            outfile.write("{}\t{}\t{}\n".format(pair[0], pair[1], score))
    outfile.close()


sentences = Sentences("testdata_pubtator_central_export.pubtator")
print(len(sentences))
word2vec = create_word2vec(sentences)
pairs = get_gene_disease_pairs("hetnet_gene_disease_pairs.tsv")
get_scores(word2vec, pairs)
print(word2vec.wv["hematologic"])


# sentences = Sentences("bioconcepts2pubtatorcentral.offset.sample")
# model = Word2Vec(sentences, size=200, window=5, min_count=1, workers=4)
# model.save("word2vec.model")
# print(model.wv.similarity("fibronectin", "vitronectin"))  # USE THIS TO SCORE EVERY GENE DISEASE PAIR
# THE ONES THAT ARE LIKELY/GROUND TRUTH SHOULD HAVE HIGHEST SCORE
# print(model.wv.doesnt_match("fibronectin vitronectine gp55 chick".split()))
# print(model.vw.most_similar(positive=['woman', 'king'], negative=['man'], topn=1))
