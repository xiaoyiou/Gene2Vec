from yxtools import dictInsert
from regulonEntity import Cis
import re

def attrCleaner(tokens):
    """
    This function simplify the attributes vector:
    1.convert strings to numbers
    2.the evidence into W or S
    3.only TFBS sequence 
    4.convert the promoter names
    """
    # 1
    tokens[3]=int(tokens[3])
    tokens[4]=int(tokens[4])
    
    # 2
    '''
    # No need for RegulonDB 9.0 since the last one is already strong/weak
    indx=len(tokens)-1
    
    if len(tokens[indx])<=0:
        tokens[indx]='W'
        return
    if tokens[indx].find('|S|')!=-1:
        tokens[indx]='S'
    else:
        tokens[indx]='W'
    '''
    # 3
    indx=len(tokens)-3
    if len(tokens[indx])<=0:
        tokens[indx]=''
        return
    seqs=re.findall("[ATCGN]+",tokens[indx])
    assert(len(seqs)==1)
    tokens[indx]=seqs[0]
    
    # 4
    temp=re.findall("[a-zA-Z]+", tokens[9])[0]
    if not temp==tokens[9]:
        print tokens[9],temp
    tokens[9]=temp
        
def getCis(path):
    
    """
    getCis takes 1 argument as the file path to 'BindingSiteSet.txt'
    and returns the result as a tulple of two dictionaries
    1.Dictionary edges: key=TF name, data=key of attrs
    2.Dictionary attrs: key=key of attrs (row #) data= tuple of detailed attributes
    3.Dicionary promoters: key=gene name, data= key of attrs
    """
    edges=dict()
    attrs=dict()
    proms=dict()
    result=(edges,attrs,proms)
    """
    cnames=[
        'TFId','TFName','TFBSID','lend'
        'rend','strand','interactionID'
        'TU','type','promoter','pos'
        'BSSeq','evidence'
    ]
    """    
    eindx=0
    with open(path,'rb') as file:
        for line in file:
            if line.startswith('#'):
                continue
            tokens=line.split('\t')
            aseq=tokens[-3]
            attrCleaner(tokens)
            if len(tokens[-3])<=0:
                continue
            attrs[eindx]=Cis(tokens[0],tokens[2],\
                             tokens[1],tokens[9],\
                             tokens[3],tokens[4],\
                             tokens[10],tokens[11],\
                             tokens[8],tokens[13],\
                             tokens[5],aseq)
            dictInsert(edges,tokens[1],eindx)
            dictInsert(proms,tokens[9],eindx)
            eindx+=1
        file.close()

    return result



            

    
