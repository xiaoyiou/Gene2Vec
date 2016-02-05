# This file deals with the visualization of the sequences
from Bio import motifs
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import numpy as np


def summ(patterns,degree,sU=True,oU=False):
    '''
    This function summarizes the patterns into sets
    patterns lD: the found patterns
    degree: the length of the patterns
    '''
    
    result=dict()
    
    for key in patterns.keys():
        seqs=[]
        for i in range(degree):
            seqs.append([])
        gaps=[]
        dirs=[]
        for p in patterns[key]:
            gaps.append(tuple(p.dist))
            dirs.append(tuple(p.dirs))
            for i in range(degree):
                seqs[i].append(p.seqs[i])
        if sU:
            for i in range(degree):
                seqs[i]=list(set(seqs[i]))
        if oU:
            gaps=list(set(gaps))
            dirs=list(set(dirs))
        temp=[len(x) for x in seqs]
        result[key]=(len(gaps),seqs,gaps,dirs,key)

    return result


def sortSum(summs):
    values=summs.values()
    return sorted(values,key=lambda x:x[0])

def visSeq(entry,degree,name):
    '''
    This function visulizes one single entry of the format
    (count,seqes[DEGREE],gaps,dirs)
    '''
    for i in range(degree):
        instances=[]
        for seq in entry[1][i]:
            instances.append(Seq(seq))
        m=motifs.create(instances)
        m.weblogo(name+str(i)+'.png',format='PNG',stack_width='large',unit_name='probability',resolution='300',color_scheme='color_classic')

def simpleVisSeq(seqs,name):
    instances=[]
    for seq in seqs:
        instances.append(Seq(seq))
    m=motifs.create(instances)
    m.weblogo(name+'.png',format='PNG',size='large',unit_name='probability',color_scheme='color_classic',stack_width='large',resolution='300')

        
def visGaps(entry,degree,bins=None):
    if degree<2 or degree>4:
        return
    if degree==2:
        # In this case we are creating a histogram
        data=[]
        for elem in entry[2]:
            data.append(elem[0])
        fig, ax = plt.subplots()
        ax.scatter(data,np.zeros_like(data)+0, s=50, alpha=0.5)

        ax.grid(True)
        fig.tight_layout()
        plt.show()
    if degree==3:
        x=[]
        y=[]
        for elem in entry[2]:
            x.append(elem[0])
            y.append(elem[1])

        fig, ax = plt.subplots()
        ax.scatter(x,y, s=50, alpha=0.5)
        ax.grid(True)
        fig.tight_layout()
    plt.show()
 
    
    
