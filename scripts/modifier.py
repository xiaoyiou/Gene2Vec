from yxtools import dictInsert
import random
import math

def dictToList(dct):
    """
    this function transforms a dict of the following structure
    dict[key]=[data1,data2,data3,...]
    into a list of tuples
    [(key1,data1),(key1,data2),(key1,data3),...]
    """

    result=[]
    for key in dct.keys():
        for elem in dct[key]:
            result.append((key,elem))

    return result

def listToDict(lst):
    """
    restore the dictionary from list of tuples. See ref of 
    dictToList()
    
    """
    result=dict()
    for (x,y) in lst:
        dictInsert(result,x,y)

    return result

def random_remove(lst,r):
    """
    This function retains r/len(lst) of elements 
    and randomly remove 1-r/len(lst) of elements
    r should be a float type

    """
    L=len(lst)
    assert(L>0)
    assert(r<=1 and r>=0)
    K=int(math.ceil((1-r)*L))    
    result= random.sample(lst,K)
    return result

def global_remove(dct,r):
    
    """
    apply the probability to all the edges and return a 
    new network
    This only test 1 r. If multiple r want to be tested 
    should call random_remove directly from a loop
    """

    elist=dictToList(dct)
    plist=random_remove(elist,r)
    print len(elist),len(plist)
    return listToDict(plist)



    
    
    
