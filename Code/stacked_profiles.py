import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from PIL import Image

#the format of the run command is as follows: python3 [path-to-this-code] [path-to-timeseries-datafile] [pulsar-period] [pulsar-sampling-period] [no.-of-pulses-to-plot-in-stackplot] [phaserange-min] [phaserange-max] [stackplot-offset] [phase-offset]

#loading the binary file data----------------------------------------------------

timeseries_path= sys.argv[1]
f= open(timeseries_path, "rb")
time_series= np.fromfile(f, dtype= 'float32')
f.close()

#info about the pulsar --------------------------------------------------------
phase_offset=0 #default value if the last argument isnt given
P = float(sys.argv[2]) #in ms
tau = float(sys.argv[3]) #in ms
num_stacked = int(sys.argv[4])
phase_range_min=float(sys.argv[5])
phase_range_max= float(sys.argv[6])

stackplot_offset=float(sys.argv[7]) #for adjusting the stack plot by a typical value of the main pulse in the folded profile

if len(sys.argv)==9:
	phase_offset = float(sys.argv[8])

num_cycles = int(np.floor(len(time_series)*tau/P)) #ignoring the last cycle which could be incomplete
N = int(round(P/tau)) #no. of samples in a cycle
M = int(round((phase_range_max - phase_range_min) * N)) #no. of samples in the given phase range

#calculating the index in time_series of the samples having zero phase-----------

zero_phase_indices=np.zeros(num_cycles,dtype=int)
for i in range(len(zero_phase_indices)):
	zero_phase_indices[i] = int(round(i*P/tau))

#cycles array having data for individual cycles----------------------------------

cycles = np.zeros((num_cycles,N)) 
for i in range(num_cycles):
	cycles[i,:] = time_series[zero_phase_indices[i]:zero_phase_indices[i]+N]

#accounting for phase_offset---------------------------------------------------

sample_offset= round(N*phase_offset)
if phase_offset != 0:
	for i in range(num_cycles):
		 
		temp=np.copy(cycles[i,:])
		temp[sample_offset:]= cycles[i,:][:N-sample_offset]
		temp[:sample_offset] = cycles[i,:][N-sample_offset:]
		cycles[i,:]=np.copy(temp)

#pulses array having data for individual cycles only for the phase range given----

pulses = np.zeros((num_cycles,M))
j = int(float(phase_range_min*N)) #position of the first sample in the desired phase range

for i in range(num_cycles):
	pulses[i,:] = cycles[i,j:j+M]

#saving the pulses data to binary or text file----------------------------------
#bin_or_text=int(input("enter 0 or 1 for saving the pulses data in binary or text format respectively: "))
bin_or_text=0

#making the header string in a particular layout. The value for each item is at the 36th character in each line
header = "#Path:                             "+timeseries_path+"\n#Period:                           "+str(P)+"\n#Sampling Interval:                "+str(tau)+"\n#num_stacked:                      "+str(num_stacked)+"\n#phase_range_min:                  "+str(phase_range_min)+"\n#phase_range_max:                  "+str(phase_range_max)+"\n#phase_offset:                     "+str(phase_offset)+"\n#no. of cycles:                    "+str(num_cycles)+"\n#no. of samples in the phase range:"+str(M)+"\nTotal no. of lines is 10(including this line) and the numeric values at each line start from the 36th character\n"

if bin_or_text==0:
	f= open("pulses.bin", 'wb')
	f.write(bytes(header, 'utf-8'))
	#pulses.tofile(f) #pulses data following the header string
	np.save(f,pulses)

elif bin_or_text==1:
	f= open("pulses.txt", 'w')
	f.write(header)
	np.savetxt(f,pulses, fmt='%.18e') #pulses data following the header string

f.close()

#making stack plot-------------------------------------------------------------

x = np.linspace(phase_range_min, phase_range_max,M)

fig = plt.figure(figsize = (12,10))
ax = plt.axes()
ax.get_yaxis().set_visible(False)
plt.xlabel('Phase')

b=0 #index of first pulse to show

plt.title("Stack plot from pulsar file- \n"+timeseries_path+"\nPulse no. {} to {}".format(b+1,b+num_stacked))

plt.ylim((0.8,stackplot_offset*num_stacked+1))
plt.xlim((phase_range_min-0.03,phase_range_max+0.03))

for i in reversed(range(num_stacked)): #going in the opposite order
	plt.plot(x,pulses[i+b,:]+i*stackplot_offset,color='black', lw=.5,zorder=2*num_stacked-2*i) #plotting the pulses in the opposite order
	if i!=0:
		plt.fill_between(x,pulses[i+b,:]+i*stackplot_offset,pulses[0+b,:], color='white', zorder=2*num_stacked -2*i-1) #filling the space between ith plot and the lowest plot for better looking stackplot
plt.fill_between(x,pulses[0+b,:],0.01, color='white')
plt.savefig("stacked_profiles.png")

#im= Image.open("stacked_profiles.png")
#im.show()

