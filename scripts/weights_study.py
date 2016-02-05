# This script contains the experiments for studying the weights
# for clustering.


## The first step is to visualize the distributions of features
## assuming entire database is available

def findHist(result):
    ll=[]
    for key in result.keys():
        ll.append((len(result[key]),result[key]))    
    return [x[1] for x in sorted(ll)]










import experiments as exp
import visualizer as vs
from yxtools import dictInsert

netpath='netdata.p'
(edges,ciss,proms)=exp.getPickle(netpath)
#sims=exp.getPickle('dmSimsL.p')
patterns=exp.constPatternReal(ciss,proms,degree=4)
sorted_sums=vs.sortSum(vs.summ(patterns,4,oU=True))


groups=exp.__regroup(ciss)


max_group=dict()
unique=dict()
for key in groups.keys():
    if key in unique.keys():
        if ciss[key].seq in unique[key]:
            continue
    dictInsert(unique,groups[key],ciss[key].seq)
    dictInsert(max_group,groups[key],key)


max_len=0
result=0
for key in max_group.keys():
    temp=len(max_group[key])
    if temp>max_len:
        max_len=temp
        result=key

    
