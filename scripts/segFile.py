'''
This python script is used to process the raw DNA sequence(s) in to iteratable sentences
of fixed-length or variant-length words
'''
import sys
import Queue
import itertools
from threading import  Thread
import cPickle as P

def gnomeGen(path):
    """
    Iterate over the gnome, and create a generator for large file or multiple
    files. Not yet required 
    """
    return None

def fmSegGnomeWorker(gnome,offset,wl,queue):
    i=offset
    result=[]
    while i<len(gnome):
        result.append(gnome[i:i+wl])
        i+=wl
    queue.put(result)


def fmSegGnome(path,wl=3):
    """
    fmSegGnome segments the entire gnome into non-overlapping words. It initiate multiple sliding windows throught the file
    path: the file dir
    wl: the fixed length of the 
    ===============
    Examples:
    ATCGGCATT, wl=3
    list1: ATC,GGC,ATT
    list2: TCG,GCA,TT
    list3: CGG,CAT,T
    """
    entire=''
    with open(path,'rb') as file:
        for line in file:
            entire+=line.rstrip()
    offsets=range(wl)
    q=Queue.Queue()
    threads=[]
    result=[]
    for offset in offsets:
        t=Thread(target=fmSegGnomeWorker,args=(entire,offset,wl,q))
        t.start()
        threads.append(t)

    for t in threads:
        t.join()
    for _ in range(wl):
        result+=q.get()

    return result

def fGenSents(wList,sLen=10000):
    i=0
    result=[]
    if len(wList)<sLen:
        return wList
    while i<len(wList):
        result.append(wList[i:i+sLen])
        i+=sLen
    return result

def fmSegFasta(path,wl=3):
    """
    fmSegFasta treat each fasta entry as a sentence and do the 
    same segementation
    """
    result=[]
    prev=''
    sentences=[]
    with open(path,'rb') as file:
        for line in file:
            line=line.rstrip()
            if line[0]=='>':
                if len(prev)>0:
                    result.append(prev)
                    prev=''
            else:
                prev+=line
    offsets=range(wl)
    for sent in result:
        tmp=[]
        for offset in offsets:
            i=offset
            while i <len(sent):
              tmp.append(sent[i:i+wl])
              i+=wl 
            sentences.append(tmp)
    return sentences
            


if __name__=="__main__":

    path=sys.argv[1]
    opath='seged_gnome.p'
    wl=3

    if(len(sys.argv))>2:
        wl=int(sys.argv[2])
        opath=(sys.argv[3])


    #data=fGenSents(fmSegGnome(path,wl))

    
    #P.dump(data,open(opath,'wb'))
    result=fmSegFasta(path,wl=wl)

    P.dump(result,open(opath,'wb'))


    

    
