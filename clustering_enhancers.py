# -*- coding: utf-8 -*-

####################################################################################################
#### Correlation between EP interacting pairs - Clustering of the dynamics of Enhancers 	  	 ###
####################################################################################################

from __future__ import division
import numpy as np
from scipy.cluster import vq
from matplotlib import pyplot as plt
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering as AgCl
from scipy.cluster import hierarchy


#GETTING THE MATRIX OF ENHANCERS

enhcs=[]
E_list=[]
E_data=open('matrixE.txt','r')
for line in E_data:
    l=line.strip('[').strip('\n').strip(']').strip("'").strip(".").split(",")
    E_list.append(l)
E_names=E_list[0]
TF_names=E_names[0].split('\t'); TF_names=TF_names[1:]
E_list=E_list[2:]

for i in range(len(E_list)):
    for j in range(len(E_list[i])):
	E_list[i][j]=E_list[i][j].strip('"')
	if j==0:
 	    E_list[i][0]=E_list[i][0].strip("'") #[2:]
	elif E_list[i][j]==" []" or E_list[i][j]==" [":
	    E_list[i][j]=E_list[i][j].replace(" []","0.0")
	    E_list[i][j]=E_list[i][j].replace(" [","0.0")
        else:
	    E_list[i][j]=E_list[i][j].strip("'[").strip("[").strip("]")
	    E_list[i][j]=E_list[i][j].strip("'")
	    E_list[i][j]=float(E_list[i][j]) #[2:]) #float(P_list[i][j][3:])

#print P_list[0]
#Creating a dictionary with the matrix   -> 'chr1:6394929-6397819'
d_e={}
for i in range(len(E_list)):
    d_e[E_list[i][0]]=E_list[i][1:]

# Function 1. Get the centroids and labels for the clustering
def Elabeling(j,nb_cl,n):
    data=[]
    e_names=[]
    for e in d_e:
        X=d_e[e][j:j+n]
        if  X!=[0.0,0.0,0.0,0.0] and X!=[0.0,0.0,0.0]:       
            e_names.append(e)
            Y=(X-np.mean(X))/np.std(X)
            data.append(Y) 
    k_means = cluster.KMeans(nb_cl)
    k_means.fit(data) 
    centroids=k_means.cluster_centers_  
    labels=[]      
    for i in k_means.labels_:
        labels.append(i)
    return labels, data,e_names,centroids

#Function 2, produce the clusters
def Eclustering(labels, data,nbclusters):
    clusters=[[] for row in range(nbclusters)]
    count=[[] for row in range(nbclusters)]
    for i in range(len(labels)):
        x=labels[i]
        clusters[x].append(data[i])
        count[x].append([labels[i],i])
    return clusters, count

# Function 3. Integrate both function 1 and 2 in order to get the results    
def Eget_results(j,nb_cl,t,n):
    labels, arr,p_names,centr=Elabeling(j,nb_cl,n)
    cl,count=Eclustering(labels,arr,nb_cl)    
    return cl,centr
    
totalc=[]
nb_cl=5 #number of clusters desired (in this case : 5)
tf=0	
TF=TF_names[tf][:-2] #number of the TF
t2=[0,30,60,120] 	#time (in mins)
Ecl,Ecentr=Eget_results(tf,nb_cl,t2,4)    	#Get the clusters (Ecl) and the centroids (Ecentr)
Etotalc=[Ecentr]

legend=[] 

#Plot the cluster trends (centroids)
plt.figure()
colour=['b','g','r','y','k','m']
for j in range(len(Etotalc)):
    legend.append('Cluster: '+str(j+5))
    for i in range(len(Etotalc[j])):
        plt.plot(t2,Etotalc[j][i],color=colour[i])
title=str(TF)+' Clustering Trends'
#plt.legend(legend)
plt.title('Enhancers:' + str(title))
plt.xticks([0,30,60,120],['0','30','60','120'])
plt.xlabel('Time (min)')
plt.show()    
 
#Plot all the clusters, with the trends for each enhancer in each cluster

for nb in range(len(Ecl)):		#Plot different figures for each cluster
    plt.figure()
    plt.xlabel('Time (min)')
    plt.ylabel('Expression')
    for i in range(len(Ecl[nb])):
        plt.plot(t2,Ecl[nb][i])  
    plt.plot(t2,Ecentr[nb],marker='*',markersize=20,color='black')    
    plt.title('Enhancers: ' +str(TF)+ ' cluster trends')
