import numpy as np
import math

# Whats given to us: 
# a numpy array named time_series containing all the sample values
# the sampling interavl tau
# the period of the pulsar P

N = 200
tau=.7
P=4.3
time_series=np.random.randint(5,size=(N))

N=np.size(time_series)
num_bins = int(P/tau)
bin_interval = P/num_bins   				 #time span of each bin

bins = [[] for i in range(num_bins)]			#this array will contains the bins and 						        their contents. Each row is represented 								by one bin	

folded= np.zeros(num_bins)				#this array will contain the average 								value for each bin

bin_count=np.zeros(num_bins,dtype=int)  		 #this keep track of how many samples we 								  have added to each bin

	
#For every sample, we have to decide which bin it goes to. .

for i in range(len(time_series)):
	bin_index = round(((i*tau)%P)/bin_interval)   #this line decides which bin gets the 							  ith sample
	
	if bin_index == num_bins: bin_index=0		
	bins[bin_index].append(time_series[i])
	bin_count[bin_index]+=1

#-------------------------------------------------------------------------
#### THE COMMENTED CODE BELOW IS A LESS SIMPLIFIED VERSION I THOUGHT. THE IDEA IS THAT THE BIN 
#WHICH CONTAINS MORE THAN HALF PORTION OF A SAMPLE WILL BE THE BIN WHERE IT GOES
	#bin_index = math.floor(((i*tau)%P)/bin_interval)	#we have to decide whether the sample goes to bin with index bin_index or the previous bin 
	#f = (((i*tau)%P) - bin_index * bin_interval)/tau	#calculates the fraction of the sample contained in the bin with index bin_index. f could be greater than 1 in case the sample lies completely inside the mentioned bin.

	#if f >=.5:		#if more than half of the sample lies in the bin with index bin_index, then it goes there.
		#bins[bin_index,bin_count[bin_index]] =time_series[i]
		#bin_count[bin_index]+=1
		
	#else:		#else the sample should go to the bin whose index is bin_index-1
		#bins[bin_index-1,bin_count[bin_index-1]] =time_series[i]
		#bin_count[bin_index-1]+=1
	
#if bin_index is 0 and f<.5 then the sample should go to the last bin. This is taken care of here since in python the index -1 means the last element of the array. 
#---------------------------------------------------------------------------------------

bins_arr=np.array(bins)
for i in range(num_bins):
	folded[i]=np.sum(bins_arr[i])/bin_count[i]	#sum all values of each bin and divide by number of samples added to each
	
print(bins_arr)
print(bin_count)
print(folded)
