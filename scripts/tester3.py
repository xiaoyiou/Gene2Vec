import experiments as exp
import random as rd
from netInfer import averageSim,maxSim
paths=['swSims.p','swSimsL.p','dmSims.p','dmSimsL.p',\
       'lengthSims.p']
names=['swAlign','swAlign-Length','directMatch',\
       'directMatchL','lengthOnly']

### random similarity measures
'''
rdSims=dict.fromkeys(range(2952))
for i in range(2952):
    rdSims[i]=dict()
    for j in range(2952):
        rdSims[i][j]=rd.random()


'''
netpath='netdata.p'
(edges,ciss,proms,sims)=exp.loadData(netpath,paths,names)
'''
### END
ratio=0.3
names.append('random')
sims['random']=rdSims
'''
'''    
(precs,recalls,ss)=exp.compareMeasureNG2(edges,sims,ratio,names,proms,wName='EntropyWeights',ciss=ciss,link_func=maxSim)
exp.visPR(precs,recalls,names+['EntropyWeights'])
'''

### From here I want to test EntropyWeightedInference Algorithm
# The similarity matrix cannot be pre-calculated
# It is dependent on the training data


