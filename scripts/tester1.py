import regulonreader as rrd
import cPickle as pickle


path='../BindingSiteSet.txt'
pp1='netdata.p'


(edges,ciss,proms)=rrd.getCis(path)

with open(pp1,'wb') as pfile:
    pickle.dump((edges,ciss,proms),pfile)
    pfile.close()
