import experiments as exp
import random as rd
from similarity import patternSim

netpath='netdata.p'
(edges,ciss,proms)=exp.getPickle(netpath)


patterns=exp.constPattern(ciss,proms)

p_data=dict()
iid=0

for key in patterns.keys():
    for p in patterns[key]:
        p_data[iid]=p
        iid+=1



'''

p_sims=dict()
for x in p_data.keys():
    p_sims[x]=dict()
    for y in p_data.keys():
        print x,y
        p_sims[x][y]=patternSim(p_data[x],p_data[y])

'''
