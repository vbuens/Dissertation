from __future__ import division

# Get the list of topological associating domains
input=open('Domainlist.txt','r')
domain=[]
for line in input:
	domain.append(line.split('\t'))
input.close()
names=domain[0]
domain=domain[1:]	#Remove the first line of the file

# Get the binding patterns of each TF (in this case, Ahr)
transcr='Ahr'
file=str(transcr)+'.txt'
input=open(file,'r')

TF=[]
for line in input:
	TF.append(line.strip('\n').split('\t'))
input.close()

#Function to detect overlaps!
def getOverlap(a,b):
    return max(0,min(a[1],b[1])-max(a[0],b[0]))

dict={}  #Dictionary with the domains, where key: chr:st-end and value will be the TFs that bind to that domain (counts!)

#For each domain, see if there is a TFBS in it
for d in range(len(domain)):
	perc=(d*100)/len(domain)
	if d%100==0: print perc 	# Print percentage to see how long it takes
	name=str(domain[d][0])+':'+str(domain[d][1])+'-'+str(domain[d][2])
	dict[name]=[]
	Dchr=domain[d][0]; Dst=int(domain[d][1]); Dend=int(domain[d][2])
	for t in range(len(TF)):
		Tchr=TF[t][0]; Tst=int(TF[t][1]); Tend=int(TF[t][2])
		if Dchr==Tchr:
			if Dst<=Tst<=Dend and Dst<=Tend<=Dend:
				dict[name].append(TF[t][3])  	#Save the TFBS in one TAD

print len(dict)

# Save results into a file
file='results/'+str(transcr)+'txt'
output=open(file,'w')				

for i in range(len(domain)):
	name=str(domain[i][0])+':'+str(domain[i][1])+'-'+str(domain[i][2])
	counts=dict[name]
	bind='\t'						#Print the results in a special format to facilitate the posterior work
	for j in range(len(counts)):
		bind=bind+str(counts[j])+'\t'
	output.write('DOMAIN:\t'+ str(name) + '\t'+'TFs:'+str(bind)+'\n')
output.close()
