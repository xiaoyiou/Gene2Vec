import math
from scipy.stats import entropy

class Pattern:
    '''
    Pattern is a tuple of Cis elements with positions,
    directions and inter-distance
    '''
    degree=0 # Number of cis elements
    ciss=None # Only store the id of the cis elements
    seqs=None  # Only store sequences
    cpos=None # The left,right positions in pairs
    dirs=None # The directions of the cis elements
    dist=None # the inter-distances, should have degree-1 items
    promName='' # promName is the promoter name of the pattern
                # patterns on the same pattern shouldn't be compared

    def __init__(self,ciss,promName):
        '''
        Parameters
        -------------
        ciss (OS): object set with the cis element 
                   object in order
        '''
        prev=None
        self.seqs=[]
        self.cpos=[]
        self.dirs=[]
        self.dist=[]
        self.ciss=[]
        for cis in ciss:
            self.degree+=1
            self.seqs.append(cis.seq)
            self.cpos.append((cis.lend,cis.rend))
            self.dirs.append(cis.strand)
            if not prev==None:
                self.dist.append(cis.lend-prev.rend)
            prev=cis
            self.ciss.append(cis)
        self.promName=promName

    def __str__(self):
        result=''
        for i in range(self.degree):
            result+='%s (%s) ' % (self.seqs[i],self.dirs[i])
            if i <self.degree-1:
                result+=str(self.dist[i])+' '
        return result

class PWM:
    '''
    PWM is a matrix that records the frequency of each base
    After initialization, the fixed-length sequences should be added
    getWeights should only be called after all sequences are added
    '''
    length=0
    freqM=None
    seqNum=0
    seqs=None
    def __init__(self, length=0,freqM=None,seqNum=0,seqs=None):
        self.length=length
        self.freqM=None
        self.seqNum=0
        self.seqs=seqs
        
    def addSeq(self,seq):
        if self.length==0:
            self.length=len(seq)
        if self.seqNum==0:
            self.seqs=[]
        self.seqs.append(seq)
        self.seqNum+=1

    def __calMatrix(self):
        assert(self.length>0 and self.seqNum>0 \
               and not self.seqs==None)
        self.freqM=[{'A':0,'T':0,'C':0,'G':0} for k in range(self.length)]
        for seq in self.seqs:
            tmpL=len(seq)
            if not tmpL==self.length:
                print "error: unequal length of seqs"
                return None
            for i in range(self.length):
                self.freqM[i][seq[i]]+=1
       

    def getPWM(self):
        if self.freqM==None:
            self.__calMatrix()
        result=[{'A':0,'T':0,'C':0,'G':0} for k in range(self.length)]
        for i in range(self.length):
            for k in self.freqM[0].keys():
                result[i][k]=self.freqM[i][k]/float(self.seqNum)
        return result
                
    def getEntropyWeights(self):
        if self.freqM==None:
            self.__calMatrix()
        upper=math.log(4)
        # 4 different base pairs
        return map(lambda x:upper-entropy(x.values()),self.freqM)
        
        


        
    
class Cis:
    'Cis is a class for cis elements'
    
    tfid=None
    bsid=None
    tfname=''
    promname=''
    lend=0
    rend=0
    pos=0.0
    seq=''
    rtype=''
    evidence=''
    strand=''
    slen=0
    def __init__(self,tfid,bsid,tfname,promname,lend,\
                 rend,pos,seq,rtype,evidence,strand,aseq=""):
        self.tfid=tfid
        self.bsid=bsid
        self.tfname=tfname
        self.promname=promname
        self.lend=lend
        self.rend=rend
        self.pos=pos
        self.seq=seq
        self.rtype=rtype
        self.evidence=evidence
        self.strand=strand
        self.slen=len(seq)
        # Adding the proximal sequences also
        self.aseq=aseq
        
    def getLen(self):
        ''' get the length of the sequence'''
        return self.slen

    def __str__(self):
        return "%s => %d %s (%d) %d %s %s %s" %\
            (self.tfname,self.lend,self.seq,\
             self.slen,self.rend,self.rtype,self.strand,self.evidence)
 
        

class BS:
    '''BS (Binding Site) is a collection of Cis elements '''
    ciss=None
    tfid=None
    tfname=''
    pwm=None
    
    def __init__(self,tfid,tfname):
        '''
        Arguments:
        tfid is the transcription factor ID assigned by RegulonDB
        tfname is gene name of tf with first letter as uppercase
        blen is the length of all the cis elements default=0
        ciss is the collection of all cis elements default=Empty
        '''
        self.tfid=tfid
        self.tfname=tfname
        self.cisNum=0
        self.pwm=PWM()
        
        
    def addCis(self,cis):
        if self.ciss==None:
            self.ciss=[]
        self.ciss.append(cis)
        self.cisNum+=1
        self.pwm.addSeq(cis.seq)


    def getPWM(self):
        return self.pwm.getPWM()

    def getEntropyWeights(self):
        return self.pwm.getEntropyWeights()

    

    
    
        
    
    
