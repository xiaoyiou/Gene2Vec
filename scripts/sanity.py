import numpy as np
import matplotlib.pyplot as plt
import math


# seperate similarity matrix
def sepSims(ciss,M):
    L=len(M)
    same_cis=[]
    dif_cis=[]
    for i in range(L):
        for j in range(i+1,L):
            if ciss[i].tfid==ciss[j].tfid and \
               ciss[i].slen==ciss[j].slen:
                same_cis.append(M[i][j])
            else:
                dif_cis.append(M[i][j])

    return (same_cis,dif_cis)

