# This file examines the relationships among the cis elements
from similarity import dmSim


def regionQuery(P,D,eps,sims):

    '''
    Parameters
    -------------
    P: the center point for the region (index)
    D: the set of all nodes (index)
    eps: the threshold of the similarity
    sims: the similarity matrix between data points

    Output
    ------------
    returns the set of P's eps-neighborhood
    '''

    result=[P]
    for cis in D:
        if sims[P][cis]>=eps:
            result.append(cis)

    return result

def regionQuery2(P,D,eps,sims):

    '''
    Parameters
    -------------
    P: the center point for the region (index)
    D: the set of all nodes (index)
    eps: the threshold of the similarity
    sims: the similarity matrix between data points

    Output
    ------------
    returns the set of P's eps-neighborhood
    '''

    result=[P]
    for cis in D:
        if sims[P][cis]>=eps:
            result.append(cis)

    return result



def expCluster(P,Nb,D,C,Ci,eps,minPts,visited,sims):
    '''
    Parameters
    ----------
    P(iN): The center point for the region (index)
    D(iS): the set of all nodes (indexes)
    C(iD): the cluster
    Ci (iN): current cluster number
    eps(rN): the similartity threshold
    Nb(iS): the neighborhood of P
    minPts(rN): the minimum points threshold
    visited(iS): the record of visted points 
    sims(DD): the similarity matrix
    '''

    for cis in Nb:
        if not cis in visited:
            visited.append(cis)
            nb=regionQuery2(cis,D,eps,sims)
            if len(nb)>=minPts:
                Nb=list(set(Nb+nb))
        if not cis in C.keys():
            C[cis]=Ci

def DBSCAN(D,eps,minPts,sims):
    visited=[]
    C=dict()
    Ci=0
    noise=[]
    for cis in D:
        if cis in visited:
            continue
        visited.append(cis)
        Nb=regionQuery2(cis,D,eps,sims)
        if len(Nb)<minPts:
            noise.append(cis)
        else:
            expCluster(cis,Nb,D,C,Ci,eps,minPts,visited,sims)
            Ci+=1

    return (noise,C)

            
    
    
