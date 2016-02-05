import cPickle as pickle
from similarity import dmSim
pfile=open('netdata.p','rb')
(edges,ciss,proms)=pickle.load(pfile)
pfile.close()

sims=dict()

for x in ciss.keys():
    sims[x]=dict()
    for y in ciss.keys():
        sims[x][y]=dmSim(ciss[x].seq,ciss[y].seq)


pfile=open('dmSims.p','wb')
pickle.dump(sims,pfile)
pfile.close()

