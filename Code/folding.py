import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from PIL import Image


#loading the binary file data----------------------------------------------------

timeseries_path= sys.argv[1]
f= open(timeseries_path, "rb")
time_series= np.fromfile(f, dtype= 'float32')
f.close()

#display_window = 20000 #no. of samples
#plt.plot(time_series[:display_window])
#plt.savefig("time_series.png")
#plt.close()	

#Folding the time series---------------------------------------------------------

#info needed for folding
phase_offset=-1 #default value if 3rd argument isnt given
P = float(sys.argv[2]) #in ms
tau = float(sys.argv[3]) #in ms
if len(sys.argv)==5:
	phase_offset = float(sys.argv[4])

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

bins_arr=np.array(bins, dtype=object)
for i in range(num_bins):
	folded[i]=np.sum(bins_arr[i])/bin_count[i]	#sum all values of each bin and divide by number of samples added to each

#print(folded)
if phase_offset != -1: #adjusting for the phase_offset in the folded profile
	bin_offset= round(num_bins*phase_offset)
	temp=np.copy(folded)
	temp[bin_offset:]= folded[:num_bins-bin_offset]
	temp[:bin_offset] = folded[num_bins-bin_offset:]
	folded=np.copy(temp)
	
#print(len(folded))
plt.plot(folded, marker='|')
plt.savefig("folded.png")
f= open("folded.txt", 'w')
np.savetxt(f,folded, fmt='%.18e')

im= Image.open("folded.png")
im.show()
