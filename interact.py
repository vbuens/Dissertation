from __future__ import division

input=open('number_of_interactions.txt','r')
interact=[]
for line in input:
	interact.append(line.strip('\n').split('\t'))
input.close()
names=interact[0]
interact=interact[1:]  	#Remove the first line of the file


# Get the binding patterns of each TF (in this case, Ahr)
transcr='Ahr'
file=str(transcr)+'.txt'
input=open(file,'r')

TF=[]
for line in input:
	TF.append(line.strip('\n').split('\t'))
input.close()


#Function to detect overlaps
def getOverlap(a,b):
    return max(0,min(a[1],b[1])-max(a[0],b[0]))

dict={}  #Dictionary with the interacts, where key: chr:st-end and value will be the TFs that bind to that interact (counts!)

#For each region, see if there is a TFBS in it
for d in range(len(interact)):
	perc=(d*100)/len(interact)
	if d%100==0: print perc  	# Print percentage to see how long it takes
	name=interact[d][3]
	dict[name]=[]
	Dchr=interact[d][0]; Dst=int(interact[d][1]); Dend=int(interact[d][2])
	for t in range(len(TF)):
		Tchr=TF[t][0]; Tst=int(TF[t][1]); Tend=int(TF[t][2])
		if Dchr==Tchr:
			if Dst<=Tst<=Dend and Dst<=Tend<=Dend:
				dict[name].append([TF[t][3],interact[d][4]])  #Save the TFBS that are located in a given region

				
#Print results in a output file
file='results/'+str(transcr)+'.txt'

output=open(file,'w')				

for i in range(len(interact)):
	name=str(interact[i][3]) 	
	counts=dict[name]
	bind='\t'					 #Save the results in a special format to facilitate the posterior work
	for j in range(len(counts)):
		bind=bind+str(counts[j][0])+','+str(counts[j][1])+'\t'
	output.write('interact:\t'+ str(name) + '\t'+'TFs:'+str(bind)+'\n')
output.close()
