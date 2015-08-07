# -*- coding: utf-8 -*-

################################################################################
#### Obtaining number of TFBSs per TAD - Clustering of TFBSs along TADs 	 ###
################################################################################

from __future__ import division
import os
import pylab as P
import matplotlib.pyplot as plt
import numpy as np

TFnames=['Ahr','Atf3','Cebpb','Ctcf','E2f1','E2f4','Egr1','Egr2','Ets','Hif1a','Irf1','Irf2','Irf4','Junb','K27Ac','Maff','Nfkb','PolII','PU1','Rel','Rela','Relb','Stat1','Stat2','Stat3']
tfn=0 	#Go through all the TFs
transcr=TFnames[tfn]  ; tfn=tfn+1
filee=str(transcr)+'.txt' 	#Open the file with the binding sites of each TFs in each TAD
inputt=open(filee,'r')
data=[]
#Read the file
for line in inputt: 		
    data.append(line.strip('\n').split('\t'))
inputt.close()
    
#Get the number of TFs per TAD
BSperTAD=[]     				
for i in range(len(data)):
    data[i]=data[i][:-1]
    tfs=data[i][3:]
    BSperTAD.append(len(tfs)) 		#append the number of TFBSs per TAD
    
#Get the number of TFs per TAD (RANDOM files previously created):
filee='Random/Random_'+str(transcr)+'.txt'
inputt=open(filee,'r')
Rdata=[] #Random
for line in inputt:
    Rdata.append(line.strip('\n').split('\t'))
inputt.close()
    
RBSperTAD=[]
for i in range(len(Rdata)):
    Rdata[i]=Rdata[i][:-1]
    tfs=Rdata[i][3:]
    RBSperTAD.append(len(tfs))    #append the number of TFs per TAD

#Obtain graphs
plt.figure()
n, bins, patches = P.hist(BSperTAD, bins=5, normed=1, histtype='stepfilled',color='blue')
Rn, Rbins, Rpatches = P.hist(RBSperTAD, bins=5, normed=1, histtype='stepfilled',color='green')
plt.clf() ; width = bins[1] - bins[0]  ;  Rwidth = Rbins[1] - Rbins[0]
plt.plot(bins[:-1],n,color='blue') 
plt.plot(Rbins[:-1],Rn,color='green') ; plt.ylim(0)
plt.plot([np.mean(BSperTAD)]*(len(bins)-1),n-0.001,color='b',linestyle='dashed')
plt.plot([np.mean(RBSperTAD)]*(len(Rbins)-1),Rn-0.001,color='g',linestyle='dashed')
plt.xlabel('Number of binding events per TADs') ; plt.ylabel('')
plt.legend([str(transcr),'Random','Mean-'+str(transcr),'Mean-Random'])
plt.title('Distribution of '+str(transcr)+' binding events. p-value: 0.9995 ')# r:'+str(round(p[tfn],4)))
title='RESULTS/'+str(transcr)+'_BSperTAD_p.jpg'
plt.savefig(title)
plt.show()    

#####################################################################
####                GET STATISTICS OF TF FILES                   ####
#####################################################################

print transcr
print min(BSperTAD)
print max(BSperTAD)
print np.mean(BSperTAD)
print np.median(BSperTAD)
    