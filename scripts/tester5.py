import modifier as mdf
from dmTools import dictThreshPREval as prEval
from yxtools import flatenDict, flatenDictList
from netInfer import netAlignInfer,getConfusion,averageSim,maxSim
import numpy as np


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
def getRecall(gold,scores,filt,thresh):
    
    filt_gold=set(gold.keys()).intersection(set(filt))
    pedges=set([key for key in scores.keys() if scores[key]>=thresh])
    com=pedges.intersection(filt_gold)
    return len(com)/float(len(filt))



# Repeats of experiment
N=1

threshs=np.arange(0,1,0.3)
precs=[]
recalls=[]

for j in range(3):
    precs.append({key:0 for key in range(len(threshs))})
    recalls.append({key:0 for key in range(len(threshs))})

ee=flatenDictList(edges)

# Create the filter list of the hidden edges
# so that the evaluation only focus on these edges
for exp_i in range(N):
    e=mdf.global_remove(edges,0.3)
    # hid_e is a list of real edges hidden
    entire_keys=edges.keys()
    subset_keys=e.keys()
    tmp_e=flatenDictList(e)
    hid_e=set(ee.keys())-set(tmp_e.keys())
    M=len(hid_e) # number of the hidden edges
    print "start inferring"
    scores1=netAlignInfer(maxSim,e,nwSims)
    scores2=netAlignInfer(maxSim,e,nwSimsL)
    scores3=netAlignInfer(maxSim,e,dmSimsL)
    ss1=flatenDict(scores1)
    ss2=flatenDict(scores2)
    ss3=flatenDict(scores3)
    for i in range(len(threshs)):
        print "thresh is ", threshs[i]
        thresh=threshs[i]
        (precs[0][i],ph)=prEval(ee,ss1,thresh)
        (precs[1][i],ph)=prEval(ee,ss2,thresh)
        (precs[2][i],ph)=prEval(ee,ss3,thresh)
        recalls[0][i]=getRecall(ee,ss1,hid_e,thresh)
        recalls[1][i]=getRecall(ee,ss2,hid_e,thresh)
        recalls[2][i]=getRecall(ee,ss3,hid_e,thresh)
    print "end of experiement", exp_i
        
    
