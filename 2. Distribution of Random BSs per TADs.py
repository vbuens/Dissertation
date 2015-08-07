# -*- coding: utf-8 -*-

##################################################################################################################
#### Create Files of random distributions of binding sites of each TF per TAD - Clustering of TFBSs along TADs ###
##################################################################################################################

from __future__ import division
import random

inputt=open('Domainlist.txt','r')
domain=[]
for line in inputt:
	domain.append(line.strip('\n').split('\t'))
inputt.close()
domain=domain[1:] 	#Remove the first line

TFnames=['Ahr','Atf3','Cebpb','Ctcf','E2f1','E2f4','Egr1','Egr2','Hif1a','Irf1','Irf2','Irf4','Junb','K27Ac','Maff','Nfkb','PolII','PU1','Rel','Rela','Relb','Stat1','Stat2','Stat3']
tfn=0 		#Go through all the TFs
transcr=TFnames[tfn]  ; tfn=tfn+1
filee=str(transcr)+'.txt'
inputt=open(filee,'r')
TF=[]
for line in inputt: 	
    l=line.strip('\n').split('\t')
    TF.append(l)  		#Get all the binding sites for a given TF
inputt.close()

dictt={}  #Dictionary with the domains, where key: chr:st-end and value will be the TFs that bind to that domain (counts!)
for d in range(len(domain)):
    name=str(domain[d][0])+':'+str(domain[d][1])+'-'+str(domain[d][2])
    dictt[name]=[]    
	
#For each TFBS, add it to a random TAD
for t in range(len(TF)):
    nb=random.randint(0,len(domain)-1) 	#Select randomly a TAD
    name=str(domain[nb][0])+':'+str(domain[nb][1])+'-'+str(domain[nb][2])
    dictt[name].append(TF[t])	 	#Assign that binding sites to a TAD
    
# Print the results into a file, similar to the one obtained of the TFBS per TAD

file='Random/Random_'+str(transcr)+'.txt'    
output=open(file,'w')				

for i in range(len(domain)):
	name=str(domain[i][0])+':'+str(domain[i][1])+'-'+str(domain[i][2])
	counts=dictt[name]
	bind='\t'
	for j in range(len(counts)):
		bind=bind+str(counts[j])+'\t'
	output.write('DOMAIN:\t'+ str(name) + '\t'+'TFs:'+str(bind)+'\n')
output.close()

#####################################################################
####                STATISTICS OF RANDOM FILES                   ####
#####################################################################

#Get the statistics of the random files just created

transcr='Ahr'
filee='Random/Random_'+str(transcr)+'.txt'
inputt=open(filee,'r')
data=[]
#Get file and store it into a variable called data
for line in inputt:
    data.append(line.strip('\n').split('\t'))
inputt.close()


RTFperTAD=[]
for i in range(len(data)):
    data[i]=data[i][:-1]
    tfs=data[i][3:]
    RTFperTAD.append(len(tfs))  #Nb of TFBSs per TAD

# Print Min, Max, Mean and Median of binding sites per TAD
print transcr
print min(RTFperTAD)
print max(RTFperTAD)
print np.mean(RTFperTAD)
print np.median(RTFperTAD)
