import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from multiprocessing import shared_memory

#loading the time series from the shared memory------------------------------------

N=1000000

shm=shared_memory.SharedMemory(name='shm_folding')
time_series=np.ndarray((N,),dtype=np.float32,buffer=shm.buf)

#Folding the time series-----------------------------------------------------------

#info needed for folding
P = 500.12345 #in ms
tau = 0.234 #in ms

num_bins = int(P/tau)					#native resolution	
bin_interval = P/num_bins   				 #time span of each bin

bins = [[] for i in range(num_bins)]			 #this array will contains the bins and their contents. Each row is represented by one bin

folded= np.zeros(num_bins)				#this array will contain the average value for each bin

bin_count=np.zeros(num_bins,dtype=int)  		 #this keep track of how many samples we have added to each bin

	
#For every sample, we have to decide which bin it goes to.
for i in range(len(time_series)):
	bin_index = round(((i*tau)%P)/bin_interval)   #this line decides which bin gets the ith sample
	
	if bin_index == num_bins: bin_index=0		
	bins[bin_index].append(time_series[i])
	bin_count[bin_index]+=1

bins_arr=np.array(bins)
for i in range(num_bins):
	folded[i]=np.sum(bins_arr[i])/bin_count[i]	#sum all values of each bin and divide by number of samples added to each

shm.close()
shm.unlink()

plt.plot(folded)
plt.show()
