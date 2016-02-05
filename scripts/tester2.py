#This file reads the network data and stores different scores
import cPickle as pickle
from similarity import swSim,nwSim,dmSim

mscore=2
misscore=-1

pp1='swSims.p'
pp2='swSimsL.p'
pp3='dmSimsL.p'



sims1=dict()
sims2=dict()
sims3=dict()


pfile_lst=[pp1,pp2,pp3]
data_lst=[sims1,sims2,sims3]

pfile=open('netdata.p','rb')
(edges,ciss,proms)=pickle.load(pfile)
pfile.close()



for x in ciss.keys():
    sims1[x]=dict()
    sims2[x]=dict()
    sims3[x]=dict()
    for y in ciss.keys():
        sims1[x][y]=swSim(ciss[x].seq,ciss[y].seq,mismatch=misscore,\
                          match=mscore)
        sims2[x][y]=swSim(ciss[x].seq,ciss[y].seq,mismatch=misscore,\
                          match=mscore,useLength=True)
        sims3[x][y]=dmSim(ciss[x].seq,ciss[y].seq,useLength=True)
for i in range(3):
    pfile=open(pfile_lst[i],'wb')
    pickle.dump(data_lst[i],pfile)
    pfile.close()



