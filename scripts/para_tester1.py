import experiments as exp
import random as rd
from similarity import patternSim
import multiprocessing
import itertools


def comp(start,end,result,data):
    flag=True
    for x in range(start,end):
        for y in range(start,end):
            if flag and x>start+500:
                print x
                flag=False
            result[(x,y)]=(patternSim(data[x],data[y]))
            


def f1(data):
    manager = multiprocessing.Manager()
    results = manager.dict()
    pool = multiprocessing.Pool(processes=8)
    N=len(data.keys())
    starts=range(0,8500,1000)
    ends=[]
    for e in starts:
        ends.append(min(e+1000,N))
        
    for (i) in range(len(starts)):
        pool.apply_async(comp, [starts[i],ends[i],results, data])
    pool.close()
    pool.join()
    return results


netpath='netdata.p'
(edges,ciss,proms)=exp.getPickle(netpath)
patterns=exp.constPattern(ciss,proms)
p_data=dict()
iid=0
for key in patterns.keys():
    for p in patterns[key]:
        p_data[iid]=p
        iid+=1


result=f1(p_data)

