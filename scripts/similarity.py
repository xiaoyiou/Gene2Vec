from skbio.alignment import StripedSmithWaterman as swAlign
from skbio.alignment import global_pairwise_align_nucleotide as nwAlign
import math
import numpy as np
from regulonEntity import PWM
from scipy.spatial import distance
def swSim(seq1,seq2,mismatch=-3,match=2,useLength=False):
    '''
    swSim uses Smith-Waterman alignment algorithm to 
    compare the sequences using local alignment 
    
    Parameters
    -----------
    seq1,seq2:  pair of sequences to be compared
    mismath(default=-3): is the mismatch score for SW alignment
    match(default=2): is the match score for SW alignment
    useLength: determines wether to ignore sequence pairs with unequal lengths
    
    Returns
    ---------
    The returned similarity value is a normalized [0,1]

    Notes
    ----------
    The gap_open_penalty is 5 and gap_extend_penalty is 2
    which are default values in skbio
    '''
    if useLength:
        if not len(seq1)==len(seq2):
            return 0

    s1=swAlign(seq1,mismatch_score=mismatch,match_score=match)
    result=s1(seq2)
    return result.optimal_alignment_score/(match*float(min(len(seq1),len(seq2))))


def nwSim(seq1,seq2,mismatch=-2,match=2,useLength=False):
    '''
    swSim uses Smith-Waterman alignment algorithm to 
    compare the sequences using local alignment 
    
    Parameters
    -----------
    seq1,seq2:  pair of sequences to be compared
    mismath(default=-2): is the mismatch score for SW alignment
    match(default=1): is the match score for SW alignment
    useLength: determines wether to ignore sequence pairs with unequal lengths
    
    Returns
    ---------
    The returned similarity value is a normalized [0,1]

    Notes
    ----------
    The gap_open_penalty is 5 and gap_extend_penalty is 2
    which are default values in skbio
    '''
    if useLength:
        if not len(seq1)==len(seq2):
            return 0
    score=nwAlign(seq1,seq2,mismatch_score=mismatch,\
                  match_score=match).score()

    return score/(match*min(len(seq1),len(seq2)))


def __dmSimEL(seq1,seq2,match,w1=None,w2=None):
    '''
    __dmSimEL(seq1,seq2,w1,w2) only compares the 
    sequences with same length
    '''
    weights=[]
    assert(len(seq1)==len(seq2))
    assert((w1==None and w2==None)or(not w1==None and not w2==None))
    if not w1==None:
        assert(len(w1)==len(w2))
        assert(len(w1)==len(seq1))
        weights=[w1[i]*w2[i] for i in range(len(w1))]
    else:
        weights=[1]*len(seq1)
    score=[weights[i]*match if seq1[i]==seq2[i] else 0 for i\
           in range(len(seq1))]
    return sum(score)/(match*float(sum(weights)))
    

def dmSim(seq1,seq2,match=2,w1=None,w2=None,useLength=False):
    '''
    dmSim(seq1,seq2,w1=None,w2=None,useLength=False):
    takes two sequences and returns the normalized similarity 
    
    Parameters
    ---------
    seq1,seq2: two sequence to be compared
    w1, w2: the weight vectors 
    useLength:determines whether to use length to filter out 
    pairs with unequal lengths of sequences
    '''
    if useLength:
        if not len(seq1)==len(seq2):
            return 0
    if len(seq1)==len(seq2):
        return __dmSimEL(seq1,seq2,match,w1,w2)

    s=''
    l=''
    lens=0
    lenl=0
    ws=[]
    wl=[]
    if len(seq1)<len(seq2):
        s=seq1
        lens=len(seq1)
        l=seq2
        lenl=len(seq2)
        ws=w1
        wl=w2
    else:
        s=seq2
        lens=len(seq2)
        l=seq1
        lenl=len(seq1)
        ws=w2
        wl=w1
    scores=[]
    for i in range(lenl-lens+1):
        tmp_seq1=l[i:i+lens]
        tmp_seq2=s
        tmp_w1=wl
        if not wl==None:
            tmp_w1=wl[i:i+lens]
        tmp_w2=ws
        scores.append(__dmSimEL(tmp_seq1,tmp_seq2,match,tmp_w1,tmp_w2))

    return max(scores)


def pwmSim(pwm1,pwm2):
    '''
    pwmSim returns the similarity measure between two PWM matrixs
    '''

    assert(pwm1.length>0 and pwm2.length>0)
    assert(pwm1.seqNum>0 and pwm2.seqNum>0)

    dimension=3*pwm1.length
    
    m1=pwm1.getPWM()
    m2=pwm2.getPWM()
    
    diff=0
    for i in range(pwm1.length):
        tmp=np.linalg.norm(\
                np.array(m1[i].values()[:-1])-\
                np.array(m2[i].values()[:-1]))**2
        diff+=tmp
        
    return math.sqrt(diff/float(dimension))


def lengthSim(seq1,seq2):
    '''
    lengthSim only compares the length between two sequences
    and return the 

    '''
    l1=len(seq1)
    l2=len(seq2)
    l=max(l1,l2)
    s=min(l1,l2)
    return 1-(l-s)/float(l)


def patternSim(p1,p2,sim_func=dmSim):
    if  p1.promName==p2.promName:
        return 0
    if not p1.degree==p2.degree:
        return 0

    seq_w=[]
    strand_w=math.log(p1.degree)
    dist_w=math.log(10)

    sqSim=0
    sSim=0
    
    
    for i in range(p1.degree):
        w=math.log(max(len(p1.seqs[i]),len(p2.seqs[i])))
        seq_sim=sim_func(p1.seqs[i],p2.seqs[i])
        sSim+=1 if p1.dirs[i]==p2.dirs[i] else 0
        sqSim+=w*seq_sim
        seq_w.append(w)

    sSim=sSim/p1.degree
    total=0
    d=distance.euclidean(p1.dist,p2.dist)
    dSim=math.exp(-d)

    result=(sSim*strand_w+sqSim+dSim*dist_w)/ \
        (strand_w+sum(seq_w)+dist_w)
    return result
    
