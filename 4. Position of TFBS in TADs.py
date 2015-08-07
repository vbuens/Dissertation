# -*- coding: utf-8 -*-

####################################################################################################
#### Obtaining the distances from TFBSs to the TAD boundaries - Position of TFBSs inside TADs  	 ###
####################################################################################################

from __future__ import division
import pylab as P
import matplotlib.pyplot as plt
import numpy as np
import os
from operator import truediv 
import math

TFnames=['Ahr','Atf3','Cebpb','Ctcf','E2f1','E2f4','Egr1','Egr2','Ets','Hif1a','Irf1','Irf2','Irf4','Junb','K27Ac','Maff','Nfkb','PolII','PU1','Rel','Rela','Relb','Stat1','Stat2','Stat3']
tfn=0  #Go through al the TFs
transcr=TFnames[tfn]  ; tfn=tfn+1  
filee=str(transcr)+'.txt'  #Get file
inputt=open(filee,'r')
data=[]
for line in inputt:
    data.append(line.strip('\n').split('\t'))
inputt.close()
for i in range(len(data)):
    data[i]=data[i][:-1]

#Initialise variables
distances=[]
dpos=[]
dneg=[]
normalp=[]
normaln=[]
 
#For each TAD:
for i in range(len(data)):
    TAD=data[i][1] ; TADchr=TAD.split(':')[0] ; TADst=int(TAD.split(':')[1].split('-')[0]) ;TADend=int(TAD.split(':')[1].split('-')[1])
	#Get centre of TAD        -Get total length of TAD
    TADc=(TADst+TADend)/2  ; TADtotal=TADend-TADst
	TFs=data[i][3:]
    for t in range(len(TFs)):  #For each one of the TF binding sites in each TAD:
        TF=TFs[t]; TFchr=TF.split(':')[0] ; TFst=TF.split(':')[1].split('-')[0]
        TFend=TF.split(':')[1].split('-')[1]; Tp=(int(TFst)+int(TFend))/2  #Get TF binding site
        if TADchr==TFchr:
            d1=(Tp-TADst)     ; dpos.append(d1)  #Get distances from the boundaries
            d2=(Tp-TADend)    ; dneg.append(d2)
            d1n=d1/TADtotal     ; normalp.append(d1n) #Get distances normalised by the length
            d2n=d2/TADtotal     ; normaln.append(d2n)
            distances.append([d1,d2])        

rawvector=dneg+dpos  	#Vector of distances from both boundaries
y=normaln+normalp  		#Vector of distances from both boundaries normalised
         
#Plot the distribution of positions of TFBSs inside topological domains:

plt.figure()
n, bins, patches = P.hist(y, bins=100, normed=0, histtype='step',color='blue')
plt.xticks([-0.5,0,0.5],['Centre','Boundary','Centre']) 		
plt.xlim(-0.5,+0.5) ; 
plt.title('Position of the '+str(transcr)+' binding site')
plt.show()   
    
## Estimate position  of TFs from the centre

#Same steps as before. Obtaning the data
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


dist=[]
norm=[]
for i in range(len(data)):
    TAD=data[i][1] ; TFs=data[i][3:] ; TADchr=TAD.split(':')[0] 
    TADst=int(TAD.split(':')[1].split('-')[0]) ;TADend=int(TAD.split(':')[1].split('-')[1])
    TADc=(TADst+TADend)/2 ; TADtotal=(TADend-TADst)/2
    for t in range(len(TFs)):
        TF=TFs[t]; TFchr=TF.split(':')[0] ; TFst=TF.split(':')[1].split('-')[0]
        TFend=TF.split(':')[1].split('-')[1]; Tp=(int(TFst)+int(TFend))/2
        if TADchr==TFchr:
            d=TADc-Tp   ;   dist.append(d) 				#Get distances from the centre
            n=(TADc-Tp)/TADtotal ; norm.append(n)		#Get distances from the centre normalised by length

#Plot the results
            
plt.figure()
n, bins, patches = P.hist(norm, bins=50, normed=1, histtype='step',color='blue')
plt.xticks([-1,0,1],['Boundary','Centre','Boundary'])
plt.xlim(-1,+1) ; 
plt.title('Position of the '+str(transcr)+' binding site')
plt.show()   

## DO SAME BUT FOR TSS: 
filee='P_TAD.txt'
inputt=open(filee,'r')
data=[]
for line in inputt:
    data.append(line.strip('\n').split('\t'))
inputt.close()
for i in range(len(data)):
    data[i]=data[i][:-1]

distances=[]
dpos=[]
dneg=[]
normalp=[]
normaln=[]
for i in range(len(data)):
    TAD=data[i][1] ; TFs=data[i][3:] ; TADchr=TAD.split(':')[0] 
    TADst=int(TAD.split(':')[1].split('-')[0]) ;TADend=int(TAD.split(':')[1].split('-')[1])
    TADc=(TADst+TADend)/2 ; TADtotal=TADend-TADst
 
    for t in range(len(TFs)):
        TF=TFs[t]; TFchr=TF.split(':')[0] ; TFst=TF.split(':')[1].split('-')[0]
        TFend=TF.split(':')[1].split('-')[1]; Tp=(int(TFst)+int(TFend))/2
        if TADchr==TFchr:
            d1=(Tp-TADst)     ; dpos.append(d1)
            d2=(Tp-TADend)    ; dneg.append(d2)
            d1n=d1/TADtotal     ; normalp.append(d1n)
            d2n=d2/TADtotal     ; normaln.append(d2n)
            distances.append([d1,d2])          

#Plot the results
y=normalp+normaln
plt.figure()
n, bins, patches = P.hist(y, bins=70, normed=1, histtype='step',color='blue')
plt.xticks([-0.5,0,0.5],['Centre','Boundary','Centre'])
plt.xlim(-0.5,+0.5) ; 
plt.title('Position of the TSSs in the TADs')
plt.show()  

##########################################################
###         DISTRIBUTION OF LENGTHS FOR THE TADs   		##
##########################################################

        
################################################################################################
###     POSITION OF TFS RESPECT TO TADs (taking into account the distribution of lengths)    ###
################################################################################################

# Get distribution of lengths of TAD      length=[>0,>100,>1000,>10000,>1Mb]
#   Array -> Nb of TADs with a certain length e.g.: [1000,500,100,42,0]

#Get the topological domains
filee='Domainlist.txt'
inputt=open(filee,'r')
data=[]
for line in inputt:
    data.append(line.strip('\n').split('\t'))
inputt.close()
data=data[1:]

#Get a vector of lengths to plot the results later (use it as a x-axis)
nb=0
Plengs=[]
while nb<=1000.001:
    Plengs.append(nb)
    nb=nb+1 
nb=-1000.001 ; Nlengs=[]
while nb<=0:
    Nlengs.append(nb)
    nb=nb+1 
lengs=Nlengs+Plengs
    
#Get the number of TADs with a length greater than a certain value (in kb)
nb=0
vectorTAD=[]
while nb<1000001:    
    nbTADs=0
    for i in range(len(data)):
        TADst=int(data[i][1]) ;TADend=int(data[i][2])
        TADpos=TADend-TADst     ;   TADneg=TADst-TADend
        if TADpos>nb:
            nbTADs=nbTADs+1     
    vectorTAD.append(nbTADs)
    nb=nb+1000 
len(vectorTAD)

vnTAD=[]
i=len(vectorTAD)-1
while i>=0:
    vnTAD.append(vectorTAD[i])
    i=i-1    
vectorTAD=vnTAD+vectorTAD    

#Names of all the TF:
names=['Ahr','Atf3','Cebpb','Ctcf','E2f1','E2f4','Egr1','Egr2','Ets','Hif1a','Irf1','Irf2','Irf4','Junb','K27Ac','Maff','Nfkb','PolII','PU1','Rel','Rela','Relb','Stat1','Stat2','Stat3']

### Function to get the distances for each TF, from its binding site to the boundaries of a TAD
def distances(transcr):
    #Read file of a given TFs and get its TFBSs per TAD
	transcr=transcr
    filee=str(transcr)+'.txt'
    inputt=open(filee,'r')
    data=[]
    for line in inputt:
        data.append(line.strip('\n').split('\t'))
    inputt.close()
    for i in range(len(data)):
        data[i]=data[i][:-1]
    dpos=[];    dneg=[]
    for i in range(len(data)):
        TAD=data[i][1] ; TADchr=TAD.split(':')[0];      TADst=int(TAD.split(':')[1].split('-')[0]) ; TADend=int(TAD.split(':')[1].split('-')[1])
		TFs=data[i][3:] 
        for t in range(len(TFs)):
            TF=TFs[t]; TFchr=TF.split(':')[0] ; TFst=TF.split(':')[1].split('-')[0]
            TFend=TF.split(':')[1].split('-')[1]; Tp=(int(TFst)+int(TFend))/2
            if TADchr==TFchr:
                d1=(Tp-TADst)/1000     ; dpos.append(d1)
                d2=(Tp-TADend)/1000    ; dneg.append(d2) 

    vectorTF=[]
	nb= 0 #Get ranges of distances 
    nbf=10 
	#Estimate the number of TFBSs located in that range
    while nb<=1000.001:
        nbTFs=0
        for i in range(len(dpos)):
            if nb<dpos[i]<nbf:
                nbTFs=nbTFs+1
        vectorTF.append(nbTFs)
        nb=nb+1  #0.10 #0.1 #000 #100
        nbf=nbf+1 #0.10 #0.1

    vnTF=[]        
    nb=-1010.001   #Get ranges of distances 
    nbf=-1000.001
	#Estimate the number of TFBSs located in that range
    while nbf<=0: 
        nbTFs=0					
        for i in range(len(dneg)):
            if nb<dneg[i]<nbf:
                nbTFs=nbTFs+1
        vnTF.append(nbTFs)
        nb=nb +1 
        nbf=nbf+1      
    vTF=vnTF+vectorTF  
    return vTF

#For each TF: call the function, get the vector of distances, divided by the vector of TAD lengths,
# 	and plot the results with the vector of lengths estimated previously
	
for i in names:
    vectorTF=distances(i)    
    totalvect=map(truediv,vectorTF,vectorTAD)
    plt.figure()
    plt.plot(lengs,totalvect, color='blue')
    plt.title('Position of TSS in TADs')
    plt.xticks([-500,0,500],['-500 Kb','Boundary','+500 Kb'])
    plt.xlim(-500,500) 
    plt.show()         
    