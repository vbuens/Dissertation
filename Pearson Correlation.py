# -*- coding: utf-8 -*-

################################################################################
#### Estimating the Linear Correlation of Enhancer and Promoter pairs	 	 ###
################################################################################

from __future__ import division
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr as corr
import random

#GETTING THE MATRIX OF PROMOTERS
proms=[]
P_list=[]
P_data=open('matrixP.txt','r')
for line in P_data:
    l=line.strip('[').strip('\n').strip(']').strip("'").strip(".").split(",")
    P_list.append(l)
P_names=P_list[0]
TF_names=P_names[0].split('\t'); TF_names=TF_names[1:]
P_list=P_list[2:]
P_bu=P_list

for i in range(len(P_list)):
    for j in range(len(P_list[i])):
        P_list[i][j]=P_list[i][j].strip('"')
        if j==0:
            P_list[i][0]=P_list[i][0].strip("'")
        elif P_list[i][j]==" []" or P_list[i][j]==" [":
            P_list[i][j]=P_list[i][j].replace(" []",'0.0')
            P_list[i][j]=P_list[i][j].replace(" [",'0.0')
            P_list[i][j]=float(P_list[i][j])
        else:
            P_list[i][j]=P_list[i][j].strip("'[").strip("[").strip("]")
            P_list[i][j]=P_list[i][j].strip("'")
            P_list[i][j]=float(P_list[i][j]) #float(P_list[i][j][3:])

 #Creating a dictionary with the matrix
d_p={}
for i in range(len(P_list)):
    d_p[P_list[i][0]]=P_list[i][1:]

#GETTING THE MATRIX OF ENHANCERS
enhcs=[]
E_list=[]
E_data=open('matrixE.txt','r')
for line in E_data:
    l=line.strip('[').strip('\n').strip(']').strip("'").strip(".").split(",")
    E_list.append(l)
E_names=E_list[0]
E_list=E_list[2:]

E_bu=E_list
#E_list=E_bu
for i in range(len(E_list)):
    for j in range(len(E_list[i])):
        E_list[i][j]=E_list[i][j].strip('"')
        if j==0:
            E_list[i][0]=E_list[i][0].strip("'")
        elif E_list[i][j]==" []" or E_list[i][j]==" [": #" '[]'" or E_list[i][j]==" '[]":
            E_list[i][j]=E_list[i][j].replace(" []",'0.0') #(" '[]'","0.0")
            E_list[i][j]=E_list[i][j].replace(" [",'0.0')
            E_list[i][j]=float(E_list[i][j])#(" '[]","0.0")
        else:
            E_list[i][j]=E_list[i][j].strip("'[").strip("[").strip("]")
            E_list[i][j]=E_list[i][j].strip("'")
            E_list[i][j]=float(E_list[i][j])

#Creating a dictionary with the matrix
d_e={}
for i in range(len(E_list)):
    d_e[E_list[i][0]]=E_list[i][1:]
    
#READ FILE WITH THE RELATION BETWEEN ENHANCER AND PROMOTERS

EP_data=open('EP_pairs.txt','r')

for line in EP_data:
    EP_list.append(line.strip('\n').split('\t'))
    
EP_pairs=[]
for i in range(len(EP_list)):	
    p=EP_list[i][3] ;    e=EP_list[i][1]  
    EP_pairs.append([p,e]) 		#Get vector of enhancer and promoter pairs

#NOW, WE ARE GOING TO COMPARE BOTH MATRIX (dictionaries) TO FIND A CORRELATION BETWEN THEM

i=0
pcorr=[]
for i in range(len(EP_pairs)):
    p=EP_pairs[i][0] 
    e=EP_pairs[i][1]
    if (p in d_p) and (e in d_e):    
        ppeaks=d_p[p];   
        epeaks=d_e[e];  
        pcorr.append([corr(ppeaks,epeaks),EP_pairs[i]]) 	#Get vector with correlation and EP pairs

highcorr=[]
lowcorr=[]        
midcorr=[]
# Get the number of EP pairs with high, mid or low correlation
for i  in range(len(pcorr)):	
    if pcorr[i][0][0]>0.8 or pcorr[i][0][0]<-0.8:
        highcorr.append(pcorr[i])
    elif pcorr[i][0][0]>0.5 or pcorr[i][0][0]<-0.5:
        midcorr.append(pcorr[i])
    else:
        lowcorr.append(pcorr[i])
        

##PLOT CORRELATION - RANDOM EP pairs

Rcorr=[]  #vector with correlations for random EP pairs
while len(Rcorr)<len(pcorr):
    i=random.randint(1,len(EP_pairs)-1)
    j=random.randint(1,len(EP_pairs)-1)
    p=EP_pairs[i][1]
    e=EP_pairs[j][0] 
    if (p in d_p) and (e in d_e):    
        ppeaks=d_p[p];  
        epeaks=d_e[e]; 
        if str(corr(ppeaks,epeaks)[0])!='nan':
            Rcorr.append([corr(ppeaks,epeaks),EP_pairs[i]])
            
 

Rhighcorr=[];   Rlowcorr=[]        ;   Rmidcorr=[]
# Get the number of Random EP pairs with high, mid or low correlation
for i  in range(len(Rcorr)):
    if Rcorr[i][0][0]>0.8  or Rcorr[i][0][0]<-0.8:
        Rhighcorr.append(Rcorr[i])
    elif Rcorr[i][0][0]>0.5 or Rcorr[i][0][0]<-0.5:
        Rmidcorr.append(Rcorr[i])
    else:
        Rlowcorr.append(Rcorr[i])

###PLOT BOTH DISTRIBUTION OF CORRELATIONS

Pc=[]
Rc=[]
for i in range(len(pcorr)):						#Plot R-square 
        Pc.append((pcorr[i][0][0]))**2)
		
for r in range(len(Rcorr)): 
        Rc.append((Rcorr[r][0][0])**2)

plt.figure() 
n1,binn1,p1=plt.hist(Pc, bins=100,color='blue',label='EP pairs',normed=1)
n,binn,p=plt.hist(Rc, bins=100,color='green',alpha=0.7,label='Random pairs',normed=1)
plt.legend()
plt.title('Distribution of correlations')
plt.xlabel('Pearson correlation (r2)')
plt.show()



#######################
## PRINT INTO A FILE ##
#######################

output=open('EP_closest.txt','w')
for i in range(len(EP_pairs)):
    output.write('Enhancer:'+'\t'+str(EP_pairs[i][1])+'\t'+'Promoter:'+'\t'+str(EP_pairs[i][0])+'\n')
    
output.close()    
    