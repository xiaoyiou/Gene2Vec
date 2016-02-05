# This module provides some common tools for
# general data mining tasks w/o information
# about the detailed
def dictThreshPREval(gold,scores,thresh):
    '''
    dictThreshPREval calculates the precision and recall 
    scores and gold should have the same number/type of keys.
    However it's not required
    
    Parameters
    ----------
    gold: {key:1/0} 1 means edges and 0 means no edges
    scores: {key:score assigned by the inferrer}
    thresh: the threshold 
    
    Retusns
    ---------
    (precision, recall)
    
    '''
    gold_set=set([key for key in gold.keys()\
                  if gold[key]==1])
    net_set=set([key for key in scores.keys()\
                 if scores[key]>=thresh])
    tpos_set=gold_set.intersection(net_set)
    tpos=float(len(tpos_set))
    ppos=float(len(net_set))
    pos=float(len(gold_set))
    prec=tpos/float(ppos) if ppos>0 else 0
    recall=tpos/float(pos) if pos>0 else 0
    return (prec, recall)

def getConfusion(gold,scores,thresh,N):
    '''
    getConfusion summarizes the confusion matrix 
    with give gold truth, scores and threshold
    
    Parameters
    ------------
    gold: the ground truth
    scores: the scores of the edges
    thresh: the cutting threshold
    N: the total number of data points
    '''
    gold_set=set([key for key in gold.keys()\
                  if gold[key]==1])
    net_set=set([key for key in scores.keys()\
                 if scores[key]>=thresh])

    tpos_set=gold_set.intersection(net_set)
    fpos_set=net_set-tpos_set
    fneg_set=gold_set-tpos_set

    TP=float(len(tpos_set))
    FP=float(len(fpos_set))
    FN=float(len(fneg_set))
    TN=N-TP-FP-FN

    return (TP/N,FP/N,FN/N,TN/N)
    

def getFilteredResult(gold,scores,filt):
    assert(len(filt)>0)
    N=len(filt)
    gkeys=gold.keys()
    skeys=scores.keys()
    fgold={key:gold[key] for key in filt if key in gkeys}
    fscore={key:scores[key] for key in filt if key in skeys}

    return (fgold,fscore,N)




