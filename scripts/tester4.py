# This file is only for the report purpose
import cPickle as pickle
import modifier as mdf
from dmTools import dictThreshPREval as prEval
from yxtools import flatenDict,flatenDictList
from netInfer import netAlignInfer,getConfusion,averageSim,maxSim
import matplotlib.pyplot as plt
import numpy as np

def getPickle(path):
    pfile=open(path,'rb')
    result=pickle.load(pfile)
    pfile.close()
    return result

'''
netp='netdata.p'
pp1='swSims.p'
pp2='swSimsL.p'
pp3='dmSimsL.p'

(edges,ciss,proms)=getPickle(netp)
nwSims=getPickle(pp1)
nwSimsL=getPickle(pp2)
dmSimsL=getPickle(pp3)
'''

# Repeats of experiments
N=5
ratio=0.3
x=[]
y=[]
x2=[]
y2=[]
x3=[]
y3=[]
ee=flatenDictList(edges)
for exp_i in range(N):
    e=mdf.global_remove(edges,0.3)
    scores1=netAlignInfer(maxSim,e,nwSims)
    scores2=netAlignInfer(maxSim,e,nwSimsL)
    scores3=netAlignInfer(maxSim,e,dmSimsL)
    ss1=flatenDict(scores1)
    ss2=flatenDict(scores2)
    ss3=flatenDict(scores3)

    for thresh in np.arange(0,1,0.1):
        result=prEval(ee,ss1,thresh)
        result2=prEval(ee,ss2,thresh)
        result3=prEval(ee,ss3,thresh)
        x.append(result[0])
        y.append(result[1])
        x2.append(result2[0])
        y2.append(result2[1])
        x3.append(result3[0])
        y3.append(result3[1])

plt.figure(1)
plt.scatter(x,y)

plt.figure(2)
plt.scatter(x2,y2)

plt.figure(3)
plt.scatter(x3,y3)

plt.show()





