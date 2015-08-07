# -*- coding: utf-8 -*-

####################################################################################################
#### Correlation between EP interacting pairs - Clustering of the dynamics of Promoters 	  	 ###
####################################################################################################

from __future__ import division
import numpy as np
from scipy.cluster import vq
from matplotlib import pyplot as plt
from sklearn.cluster import KMeans
from sklearn import cluster

#GETTING THE MATRIX OF PROMOTERS
proms=[]
P_list=[]
#P_data=open('GOOD/avmatrixP.txt','r')
P_data=open('matrixP.txt','r')
for line in P_data:
    l=line.strip('[').strip('\n').strip(']').strip("'").strip(".").split(",")
    P_list.append(l)
P_names=P_list[0]
TF_names=P_names[0].split('\t'); TF_names=TF_names[1:]

for i in range(len(P_list)):
    for j in range(len(P_list[i])):
	P_list[i][j]=P_list[i][j].strip('"')
	if j==0:
 	    P_list[i][0]=P_list[i][0].strip("'")
	elif P_list[i][j]==" []" or P_list[i][j]==" [":
	    P_list[i][j]=P_list[i][j].replace(" []","0.0")
	    P_list[i][j]=P_list[i][j].replace(" [","0.0")
        else:
	    P_list[i][j]=P_list[i][j].strip("'[").strip("[").strip("]")
	    P_list[i][j]=P_list[i][j].strip("'")
	    P_list[i][j]=float(P_list[i][j]) #float(P_list[i][j][3:])

#Creating a dictionary with the matrix   -> 'chr1:6394929-6397819'
d_p={}
for i in range(len(P_list)):
        d_p[P_list[i][0]]=P_list[i][1:]


# Function 1. Get the centroids and labels for the clustering
def Plabeling(j,nb_cl,n):
    data=[]
    p_names=[]
    for p in d_p:#, loc='upper left')
        X=d_p[p][j:j+n]
        if X!=[0.0,0.0,0.0,0.0]:           
            p_names.append(p)
            Y=(X-np.mean(X))/np.std(X)
            data.append(Y) 
    k_means = cluster.KMeans(nb_cl)
    k_means.fit(data) 
    centroids=k_means.cluster_centers_  
    labels=[]      
    for i in k_means.labels_:
        labels.append(i)
    return labels, data,p_names,centroids

#Function 2, produce the clusters
def Pclustering(labels, data,nbclusters):
    clusters=[[] for row in range(nbclusters)]
    count=[[] for row in range(nbclusters)]
    for i in range(len(labels)):
        x=labels[i]
        clusters[x].append(data[i])
        count[x].append([labels[i],i])
    return clusters, count

# Function 3. Integrate both function 1 and 2 in order to get the results        
def Pget_results(j,nb_cl,t,n):
    labels, arr,p_names,centr=Plabeling(j,nb_cl,n)
    cl,count=Pclustering(labels,arr,nb_cl)    
    return cl,centr
    
totalc=[]
nb_cl=5  #number of clusters desired (in this case : 5)
tf=8
TF=TF_names[tf][:-2] #number of the TF
t2=[0,30,60,120] 	#time (in mins)
cl,centr=Pget_results(tf,nb_cl,t2,4)        	#Get the clusters (cl) and the centroids (centr)
totalc=[centr]

legend=[] 

#Plot the cluster trends (centroids)
plt.figure()
colour=['b','g','r','y','k','m']
for j in range(len(totalc)):
    legend.append('Cluster: '+str(j+5))
    for i in range(len(totalc[j])):
        plt.plot(t2,totalc[j][i],color=colour[i])
title=str(TF)+' Clustering Trends'
plt.title('Promoters: ' + str(title))
plt.xticks([0,30,60,120],['0','30','60','120'])
plt.xlabel('Time (min)')
plt.show()    
    
#Plot all the clusters, with the trends for each enhancer in each cluster

for nb in range(len(cl)): 	#Plot different figures for each cluster
    title='Promoters: '+str(TF)+ ' clustering trends' #'_'+str(len(cl))+'cl_'+'nb'+str(nb+1)+'.png'
    plt.figure()
    plt.xlabel('Time (minutes)')
    plt.xticks([0,30,60,120],['0','30','60','120'])
    plt.ylabel('Expression')
    for i in range(len(cl[nb])):
        plt.plot(t2,cl[nb][i])  
    plt.plot(t2,centr[nb],marker='*',markersize=20,color='black')    
    plt.title(title)
