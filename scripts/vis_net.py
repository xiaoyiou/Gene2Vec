import experiments as exp
import modifier as mdf
import re

paths=['swSims.p','swSimsL.p','dmSims.p','dmSimsL.p',\
       'lengthSims.p']
names=['swAlign','swAlign-Length','directMatch',\
       'directMatchL','lengthOnly']

netpath='netdata.p'

(edges,ciss,proms)=exp.getPickle(netpath)

filt1=[]
filt2=[]
#ratio=0.99
#filt1=mdf.random_remove(edges.keys(),ratio)


'''
for cis in ciss.keys():
    entry=ciss[cis]
    target=entry.promname
    target= re.sub(r'p\d?$','',target).upper()
    source=entry.tfname.upper()
    if not (source,target) in filt1:
        filt1.append((source,target))


for elem in filt1:
    print elem[0],'R',elem[1]
'''

tss=exp.getTS2(edges,ciss,proms,None,0.1)
factors=exp.getCooccur(*tss,min_c=0.1)

for key in factors.keys():
    print key[0],factors[key],key[1]
