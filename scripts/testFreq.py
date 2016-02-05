import experiments as exp
from yxtools import dictInsert
import math
import matplotlib.pyplot as plt
import operator
(edges,ciss,proms)=exp.getPickle('netdata.p')

transactions=exp.getTS(edges,ciss,proms)

rules=exp.getAssoc(transactions,min_s=2,\
                   min_c=0.5)

# Determine the relationship between support and
# size of frequent itemset


(perfs,maks)=exp.evalCisCollection(rules,weighted=True)

'''
y=[[],[],[],[],[]]

for i in perfs.keys():
    for elem in perfs[i]:
        y[i-1].append(elem)
'''
x=[]
for key in maks.keys():
    x.append(len(maks[key][0])+1)

plt.hist(x)
plt.show()


'''
fig=plt.figure(1,figsize=(9,6))
ax=fig.add_subplot(111)
ax.set_xticklabels(['2','3','4','5','6'])
plt.ylim([-.1,1.1])
plt.xlabel('Length of frequent itemsets')
plt.ylabel('Confidence %')
bp=ax.boxplot(y)
plt.show()
   


plt.scatter(x,y)
plt.show()
'''
