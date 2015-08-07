# -*- coding: utf-8 -*-

##################################################################################################################
#### Create Files of random distributions of binding sites of each TF per TAD - Clustering of TFBSs along TADs ###
##################################################################################################################

from __future__ import division
import math
import matplotlib.pyplot as plt
import numpy as np
from operator import truediv 
import scipy as sc

### GET THE TOTAL NUMBER OF INTERACTIONS

filee='number_of_interactions.txt'   # read file
inputt=open(filee,'r')
data=[]
for line in inputt:
    data.append(line.strip('\n').split('\t'))
inputt.close()
#Get the total number of interactions per region
Tinterct=[]
for i in range(len(data)):
    Tinterct.append(int(data[i][4]))
    
## GET NUMBER OF INTERACTIONS PER REGION WITH TFBSs

TFnames=['Ahr','Atf3','Cebpb','Ctcf','E2f1','E2f4','Egr1','Egr2','Hif1a','Irf1','Irf2','Irf4','Junb','K27Ac','Maff','Nfkb','PolII','PU1','Rel','Rela','Relb','Stat1','Stat2','Stat3']
tfn=0   #Go through all the TFs
transcr=TFnames[tfn]  ; tfn=tfn+1
filee=str(transcr)+'.txt' #Read file
inputt=open(filee,'r')
data=[]
for line in inputt:
    data.append(line.strip('\n').split('\t'))
inputt.close()
for i in range(len(data)):
    data[i]=data[i][:-1]

#Get the number of interaction
nbinterct=[]
for i in range(len(data)):
    TFS=data[i][3:]
    for t in range(len(TFS)):
        TF=TFS[t].split(',')
        nbinterct.append(int(TF[1]))

#Plot the results
plt.figure()
n1,bins1,patch1=plt.hist(nbinterct,bins=np.logspace(0.1,10.0,100),label=str(transcr))#,normed=1)  
n2,bins2,patch2=plt.hist(Tinterct,bins=np.logspace(0.1,10.0,100),label='All interactions')#,normed=1)
plt.clf()
x1_frac = n1/float(sum(n1)) # Each bin divided by total number of objects.
width1 = bins1[1] - bins1[0] # Width of each bin.
x1 = np.ravel(zip(bins1[:-1], bins1[:-1]+width1))
y1 = np.ravel(zip(x1_frac,x1_frac))
x2_frac = n2/float(sum(n2)) # Each bin divided by total number of objects.
width2 = bins2[1] - bins2[0] # Width of each bin.
x2 = np.ravel(zip(bins2[:-1], bins2[:-1]+width2))
y2 = np.ravel(zip(x2_frac,x2_frac))

plt.plot(x1,y1,label=str(transcr)) 
plt.plot(x2,y2,label='All interactions')
plt.title('Distribution of interactions along regions with TFBSs.')
plt.legend(loc=2) #Change the location of the legend
plt.xlabel('Number of interactions per region (log)')
plt.ylabel('Number of regions')
plt.gca().set_xscale("log") 	#log scale
plt.show()


######################################################### 
## STATISTIC TEST : Z-test of comparision of means  	#
#########################################################
#Vector to store the results of the tests
sampletest=[]
def twoSampZ(X1, X2, mudiff, sd1, sd2, n1, n2):
    from numpy import sqrt, abs, round
    from scipy.stats import norm
    pooledSE = sqrt(sd1**2/n1 + sd2**2/n2)
    z = ((X1 - X2) - mudiff)/pooledSE
    pval = 2*(1 - norm.cdf(abs(z)))
    return round(z, 3), round(pval, 4)

#Estimate results for each factor transcription
tfn=0    
transcr=TFnames[tfn]  ; tfn=tfn+1
filee=str(transcr)+'.txt'
inputt=open(filee,'r')
data=[]
for line in inputt:
    data.append(line.strip('\n').split('\t'))
inputt.close()
for i in range(len(data)):
    data[i]=data[i][:-1]

nbinterct=[]
for i in range(len(data)):
    TFS=data[i][3:]
    for t in range(len(TFS)):
        TF=TFS[t].split(',')
        nbinterct.append(int(TF[1]))

x1=np.mean(nbinterct) 	#interactions per region (regions with TFBSs)
x2=np.mean(Tinterct)   	#interactions per region (all regions of the genome)
sd1=np.std(nbinterct) ; sd2=np.mean(Tinterct) 
n1=len(nbinterct)   ; n2=len(Tinterct)
x=twoSampZ(x1,x2,0,sd1,sd2,n1,n2)    
sampletest.append(x) #Get vector with the results for all the TFs