# This defines the methods to experiment binding site sequences using 
# vectorized words from Word2Vec 
import gensim
import regulonEntity
import scipy.spatial.distance as ssd

def seq2vec(seq,model,wl):
    """
    seq: the sequence waiting to be transformed
    model: the gensim W2V model to be used
    wl: the length of each word in the model
    """
    
    i=0
    result=0
    while i<len(seq):
        if len(seq)-i-1<wl:
            break
        result+=model[seq[i:i+wl]]
        i+=wl
    return result


def transCis(ciss,wl,model,proximal=False):
    """
    ciss: the dict database of Cis objects
    wl: the length of words
    model: the gensim W2V model
    --------------------------------------
    returns a dictionary of vectors, key is the cis index and value is the numpy.ndarray
    """
    
    result=dict()
    for key in ciss.keys():
        if proximal:
            result[key]=seq2vec(ciss[key].aseq.upper(),model,wl)
        else:
            result[key]=seq2vec(ciss[key].seq,model,wl)
    return result


    
        
def printSeqs(lst,ciss):
    for c in lst:
        print ciss[c]

def printVecs(lst,vCiss):
    for c in lst:
        print vCiss[c]


def pwAverageCosine(edges,vCiss,thresh=5):

    result=dict()

    for key in edges:
        if len(edges[key])<thresh:
            continue
        tempV=0
        count=0
        clst=edges[key]
        for x in clst:
            for y in clst:
                if x==y:
                    continue
                count+=1
                tempV+=1-ssd.cosine(vCiss[x],vCiss[y])
        tempV/=count
        result[key]=tempV
    return result

def pwAverageCosineAll(vCiss,thresh=5):
    clst=vCiss.keys()
    tempV=0
    count=0
    for x in clst:
        for y in clst:
            if x==y:
                continue
            count+=1
            tempV+=1-ssd.cosine(vCiss[x],vCiss[y])
    return tempV/count
    
    
    
