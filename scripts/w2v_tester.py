import gensim
import cPickle as P
import vCis
net_path='../data/netdata.p'
model_path='../data/6_inter_model.dat'


(edges,ciss,proms)=P.load(open(net_path,'rb'))
model=gensim.models.Word2Vec.load(model_path)

vciss=vCis.transCis(ciss,6,model,proximal=True)


