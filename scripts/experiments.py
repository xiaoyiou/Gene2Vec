# This file contains some handy tools for experimenting
import cPickle as pickle
import numpy as np
import modifier as mdf
import matplotlib.pyplot as plt
import random as rd
from dmTools import dictThreshPREval as prEval
from yxtools import flatenDict,flatenDictList,dictInsert
from netInfer import netAlignInfer,getConfusion,averageSim,maxSim
from regulonEntity import PWM,Pattern
from similarity import dmSim
import math
from pymining import itemmining,assocrules
def getPickle(path):
    pfile=open(path,'rb')
    result=pickle.load(pfile)
    pfile.close()
    return result

def loadData(netpath,paths,names):
    assert(len(paths)==len(names))
    (edges,ciss,proms)=getPickle(netpath)
    sims=dict.fromkeys(names)
    for i in range(len(paths)):
        sims[names[i]]=getPickle(paths[i])
        
    return (edges,ciss,proms,sims)




def compareMeasureNG(edges,sims,ratio,names,wName=None,ciss=None,link_func=maxSim):
    '''
    NG means this experiment is not using the grouping
    of cis-element
    '''

    precs=dict.fromkeys(names)
    recalls=dict.fromkeys(names)
    ss=dict.fromkeys(names)
    
    ee=flatenDictList(edges)#ee is true edges
    e=mdf.global_remove(edges,ratio)
    for name in names:
        ss[name]=flatenDict(netAlignInfer\
                            (link_func,e,sims[name]))
        
    for name in names:
        precs[name]=[]
        recalls[name]=[]
        for thresh in np.arange(0,1.0,0.1):
            (p,r)=prEval(ee,ss[name],thresh)
            precs[name].append(p)
            recalls[name].append(r)

    if not wName==None:
        (precs[wName],recalls[wName])=doEntropyWeights(ee,e,ciss,\
                                    link_func=link_func)
    return (precs,recalls,ss)

def visPR(precs,recalls,names):
    '''
    precs and recalls are dictionaries of the 
    precisions and recalls 
    '''
    assert(precs.keys()==recalls.keys())
    for name in names:
        lines=plt.plot(recalls[name],precs[name])
        plt.setp(lines, linewidth=2.0)
    plt.ylim([0,1])
    plt.xlim([0,1])
    plt.legend(names,loc='left right')
    plt.show()

def getEntropyWeights(e,ciss):
    '''
    getEntropyWeights calculates the similarity based on 
    partial edges e
    
    Retruns
    -----------
    weights: a dictionary of weight associated with each cis-element {2019:[0.1,0.4,0.9...]}  
    '''
    weights=dict.fromkeys(ciss.keys())
    groups=dict.fromkeys(e.keys())
    for key in ciss.keys():
        weights[key]=[1.0]*ciss[key].slen

    
    for key in groups.keys():
        all_cis=e[key]
        if len(all_cis)==0:
            continue
        groups[key]=dict()
        for cis in all_cis:
            cc=ciss[cis]
            dictInsert(groups[key],cc.slen,cis)
    for key in groups.keys():
        for l in groups[key].keys():
            pwm=PWM()
            for cisId in groups[key][l]:
                elem=ciss[cisId]
                pwm.addSeq(elem.seq)
            tmp_weights=pwm.getEntropyWeights()
            for cisId in groups[key][l]:
                weights[cisId]=tmp_weights
            
    return weights

def doEntropyWeights(ee,e,ciss,link_func=maxSim,proms=None):
    precs=[]
    recalls=[]
    
    weights=getEntropyWeights(e,ciss)
    cis_keys=ciss.keys()
    sims=dict()
    for x in cis_keys:
        sims[x]=dict()
        for y in cis_keys:
            sims[x][y]=dmSim(ciss[x].seq,ciss[y].seq\
                             ,w1=weights[x],w2=weights[y]\
                             ,useLength=True)
    ss=flatenDict(netAlignInfer(link_func,e,sims))
    for thresh in np.arange(0,1,0.1):
        
            (p,r)=prEval(ee,ss,thresh)
            precs.append(p)
            recalls.append(r)
    return (precs,recalls)

def __regroup(ciss):
    '''
    regroup returns a dict of the BS labels for each cis element
    the cis-elements with the same label are considered to be
    the same 

    Returns
    -----------
    results: a dictionary of relationships from cis # to \
    TFName+cis_length as the new key 
    '''

    result=dict()
    for key in ciss.keys():
        obj=ciss[key]
        rkey=obj.tfname+str(obj.slen)
        result[key]=rkey

    return result

def __findTs(edges,ciss,proms):
    '''
    findTs will give the transactions lists indexed by the TF name
    proms is the transformed dictionary of promoter
    
    Parameters
    -----------
    edges

    Returns
    ----------
    result: the dict with 
    '''
    result=dict()
    nkey=''
    for key in edges.keys():
        cis_lst=edges[key]
        prom_set=set()
        filt=[]
        for cis in cis_lst:
            if ciss[cis].promname in filt:
                continue
            filt.append(ciss[cis].promname)
            nkey=key+str(ciss[cis].slen)
            dictInsert(result,nkey,proms[ciss[cis].promname])
            
        
    return result


def __swapCis(proms,filt):
    '''
    replace the proms with new naming of the cis elements
    '''
    result=dict()
    for key in proms.keys():
        lst=proms[key]
        result[key]=list(set([filt[elem] for elem in lst]))
    
    return result

def getTS(edges,ciss,proms):
    '''
    getTS get the netdata and returns the dictionary
    with {tfName:[[FNR14,DnaA9...],[DnaA9,DpiA23...],[...]]}
    '''
    return __findTs(edges,ciss,__swapCis(proms,__regroup(ciss)))

def getAssoc(transactions,min_s=2,min_c=0.5):
    '''
    getAssoc will return the association rule in the following 
    format
    '''
    result=dict()
    for key in transactions.keys():
        relim_input=itemmining.get_relim_input(\
            transactions[key])
        itemset=itemmining.relim(relim_input\
                                 ,min_support=min_s)
        result[key]=assocrules.mine_assoc_rules(\
                                                itemset,min_support=min_s,min_confidence=min_c)

    return result


def evalCisCollection(rules,weighted=False):
    '''
    evalCisCollection stores the 
    (length,support,confidence) in the dictionary list
    '''
    results=dict()
    max_len=0
    maks=dict()
    for key in rules.keys():
        if len(rules[key])==0:
            continue
        max_len=0
        for elem in rules[key]:
            if len(elem[0])>1:
                continue
            if not elem[0].issuperset({key}):
                continue
            dictInsert(results,key,\
                       (len(elem[1]),elem[2],elem[3]))
            if len(elem[1])>max_len:
                maks[key]=(elem[1],elem[2],elem[3])
                max_len=len(elem[1])
            
    perfs=dict()
    for key in results.keys():
        for elem in results[key]:
            l=elem[0]
            if weighted:
                for i in range(elem[1]):
                    dictInsert(perfs,l,elem[2]) 
            else:
                dictInsert(perfs,l,elem[2])
    
    return (perfs,maks)
    
def getTS2(edges,ciss,proms,scores,thresh):
    '''
    getTS2 returns the transactions for association rules mining
    '''

    groups=dict()
    reverse=dict()
    for key in ciss.keys():
        groups[key]=key
    
    for key in edges.keys():
        for cis in edges[key]:
            obj=ciss[cis]
            rkey=obj.tfname+str(obj.slen)
            groups[cis]=rkey
            reverse[rkey]=cis
    if not scores ==None:
        for key in scores.keys():
            if scores[key]>=thresh:
                groups[key[1]]=key[0]+str(ciss[key[1]].slen)
            
            # At this point the name change is ready to use
    

    transactions=[]
    for key in proms.keys():
        cis_lst=proms[key]
        filt=[]
        for cis in cis_lst:
            if groups[cis] in filt:
                continue
            filt.append(groups[cis])
        transactions.append(filt)

    return (transactions,groups,reverse)

def getAssoc2(ts,groups,reverse,min_s=2,min_c=0.5):
    crules=dict()
    result=dict()
    relim_input=itemmining.get_relim_input(ts)
    itemset=itemmining.relim(relim_input,min_support=min_s)
    rules=assocrules.mine_assoc_rules(itemset,min_support=min_s\
                                      ,min_confidence=min_c)
    # Now calculate the best rule for each cis
    # Clean the rules
    for rule in rules:
        if len(rule[0])>1:
            continue
        else:
            if rule[0] in crules.keys():
                if (len(rule[1])+1)*rule[3]<=crules[rule[0]]:
                    continue
                crules[rule[0]]=(len(rule[1])+1)*rule[3]
    for cis in groups.keys():
        key=frozenset({groups[cis]})
        if key in crules.keys():
            result[cis]=crules[key]
        
    return result
def getCooccur(ts,groups,reverse,min_s=2,min_c=0.5):
    crules=dict()
    result=dict()
    relim_input=itemmining.get_relim_input(ts)
    itemset=itemmining.relim(relim_input,min_support=min_s)
    rules=assocrules.mine_assoc_rules(itemset,min_support=min_s\
                                      ,min_confidence=min_c)
    # Now calculate the best rule for each cis
    # Clean the rules
    for rule in rules:
        if len(rule[0])>1:
            continue
        else:
            if not rule[0] in crules.keys():
                crules[rule[0]]=dict()
                for elem in rule[1]:
                    crules[rule[0]][elem]=rule[3]
    for x in reverse.keys():
        if not frozenset({x}) in crules.keys():
            continue
        for y in reverse.keys():
            if not y in crules[frozenset({x})].keys():
                continue
            result[(x,y)]=crules[frozenset({x})][y]
            
        
    return result

def __getSubset(lst,remain):
    result=[]
    if remain==1:
        for elem in lst:
            result.append([elem])

    else:
        for i in range(len(lst)-remain):
            t=__getSubset(lst[i+1:],remain-1)
            for elem in t:
                result.append([lst[i]]+elem)
    return result
def constPattern(ciss,proms,degree=2):
    '''
    constPattern create new data points that are tuples of
    cis elements in order
    '''

    result=dict()
    for key in proms.keys():
        prom=proms[key]
        if len(prom)<degree:
            continue
        clst=[]
        for c in prom:
            clst.append(ciss[c])
        clst.sort(key=lambda x: float(x.lend))
        
        patterns=__getSubset(clst,degree)
        for elem in patterns:
            dictInsert(result,key,Pattern(elem,key))
    return result
      
def constPatternReal(ciss,proms,degree=2):
    '''
    constPatternReal creates groups of patterns based on
    their real labelings
    '''
    result=dict()
    temp=constPattern(ciss,proms,degree)
    for key in temp.keys():
        for p in temp[key]:
            name=''
            existing=[]
            for cis in p.ciss:
                cKey=cis.tfname
                if not cKey in existing:
                    existing.append(cKey)
                    name+=cKey
            if len(existing)<degree:
                continue
            dictInsert(result,name,p)
    return result
        

    



    
def compareMeasureNG2(edges,sims,ratio,names,proms,wName=None,ciss=None,link_func=maxSim):
  
    precs=dict.fromkeys(names)
    recalls=dict.fromkeys(names)
    ss=dict.fromkeys(names)
    
    ee=flatenDictList(edges)#ee is true edges

    e=mdf.global_remove(edges,ratio)
    for name in names:
        ss[name]=flatenDict(netAlignInfer\
                            (link_func,e,sims[name]))
    tss=getTS2(edges,ciss,proms,ss[name],0.1)
    factors=getAssoc2(*tss,min_c=0.9)
    for name in names:
        precs[name]=[]
        recalls[name]=[]
        for thresh in np.arange(0,1.0,0.1):
            
            for key in ss[name].keys():
                if key[1] in factors.keys():
                    if key[0]+str(ciss[key[1]].slen) == tss[1][key[1]]:
                        #print 'BEFORE',ss[name][key]
                        ss[name][key]=1
                        #print 'AFTER',ss[name][key]
                    
            (p,r)=prEval(ee,ss[name],thresh)
            precs[name].append(p)
            recalls[name].append(r)

    if not wName==None:
        (precs[wName],recalls[wName])=doEntropyWeights(ee,e,ciss,\
                                    link_func=link_func)
    return (precs,recalls,ss)



    
    
if __name__=='__main__':
    paths=['swSims.p','swSimsL.p','dmSims.p','dmSimsL.p',\
           'lengthSims.p']
    names=['swAlign','swAlign-Length','directMatch',\
           'directMatchL','lengthOnly']

    ### random similarity measures
    #'''
    rdSims=dict.fromkeys(range(2952))
    for i in range(2952):
        rdSims[i]=dict()
        for j in range(2952):
            rdSims[i][j]=rd.random()
    #'''
    ### END
    
    netpath='netdata.p'
    ratio=0.3
    (edges,ciss,proms,sims)=loadData(netpath,paths,names)
    names.append('random')
    sims['random']=rdSims

    (precs,recalls)=compareMeasureNG(edges,sims,ratio,names)
    visPR(precs,recalls,names)
    #visPR(precs,recalls,'lengthOnly')

