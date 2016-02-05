# This script is testing the ideas of the w2v experiments
# using ecoli entire gnome
import cPickle as P
import regulonreader as rrd
import gensim

sent_path='../data/4ecoli.p'
net_path='../data/netdata.p'

sentences=P.load(open(sent_path,'rb'))
(edges,ciss,proms)=P.load(open(net_path,'rb'))


#----------Training of W2V model-----

#Some parameters of the training
wd=20 # looking at total of 120 base pairs
wk=8 # number of threads
size=30 # dependent on the vocab size, we have 4^4=265
min_count=10 # not very useful beacuse genes have a smaller vocabluary, so repeats are more probable
sg=0 # CBOW or Skip-Gram
hs=0 # hierarchical softmax or not
it=15 #number of iterations (epochs) over the gnome


model=gensim.models.Word2Vec(sentences,window=wd,
 workers=wk,size=size, min_count=min_count, sg=sg, 
 hs=hs, iter=it)





