# This module provides algorithms for netInference algorithms
import numpy as np
from regulonEntity import PWM,Cis


def evaluate(gold,scores,thresh):
    """
    This function returns the confusion matrix as 
    the result
    gold should be in the following format:
    [(TF,target1),(TF,target2),...]
    and scores should be in the following format:
    {(TF,target1):score, (TF,target2):score,...}
    """
    N=len(gold)
 
    gold_set=set(gold)
 
    net_set=set([elem for elem in scores.keys()\
                 if scores[elem]>=thresh])
    retrived=len(net_set)
    relevant=len(gold_set)
    common=float(len(gold_set.intersection(net_set)))
    prec=common/retrived
    recall=common/relevant
    return (prec,recall)
    


def getConfusion(edges,scores,thresh):
    gold=[]
    scores_f=dict()
    for key in edges.keys():
        if not key in scores.keys():
            continue
        for target in edges[key]:
            gold.append((key,target))

    for key in scores.keys():
        for target in scores[key].keys():
            scores_f[(key,target)]=scores[key][target]
    return evaluate(gold,scores_f,thresh)



def maxSim(group_ids,target_id,sims):
    if len(group_ids)==0:
        return 0
    scores=[]
    for idx in group_ids:
        if not idx==target_id:
            scores.append(sims[idx][target_id])
    if len(scores)==0:
        return 0
    else:
        return max(scores)


def averageSim(group_ids,target_id,sims):
    if len(group_ids)==0:
        return 0
    scores=[]
    for idx in group_ids:
        scores.append(sims[idx][target_id])
    return sum(scores)/float(len(scores))


def netAlignInfer(scoreFunc,edges,sims):
    '''
    netAlignInfer use existing similarity matrix between pairs of cis elements
    Given that no groups of Cis are required, so we can pre-calculate the similarity
    matrix

    Parameters
    ---------
    edges: dictionary of existing edges edges[key]=[cis_id1,cis_id2...]
    sims: is a pairwise similarity matrix of Cis elements
    scoreFunc: determines similarity between C' and cis where C' is a set of Cis
    Returns
    ---------
    A dictionary of scores assigned to TF-BS relationships
    result[TF_name][cis_id] in range [0,1]

    '''
    result=dict()
    
    for key in edges.keys():
        result[key]=dict.fromkeys(sims.keys(),0)
        lst_motifs=edges[key]
        n=len(lst_motifs)
        if n<=0:            
            continue
        for target in sims.keys():
            result[key][target]=scoreFunc(lst_motifs,target,sims)
        
    return result



def netWeightedInfer(edges,ciss,s=2.0,f=-0.1):
    result=dict()
    groups=dict()#key is cis_id, value is pwm_id 
    pwms=dict()# key is pwm_id, value is PWM instance
    weights=dict() # key is cis_id, value is weight vector

    # Creating PWMs for weightCalculation
    pwm_id=0

    for key in edges.keys():
        result[key]=dict.fromkeys(sims.keys(),0)
        lst_motifs=edges[key]
        n=len(lst_motifs)
        if n<=0:            
            continue
            
        pwms[pwm_id]=PWM()
        for indx in lst_motifs:
            pwms[pwm_id].addSeq(ciss[indx].seq)
        
        
            
